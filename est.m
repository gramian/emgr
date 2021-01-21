function [r,m] = est(sys,task,config)
%%% project: emgr - EMpirical GRamian Framework ( https://gramian.de )
%%% version: 5.9 (2021-01-21)
%%% authors: Christian Himpe (0000-0003-2194-6754)
%%% license: BSD-2-Clause (opensource.org/licenses/BSD-2-Clause)
%%% summary: est - empirical system theory (emgr frontend)

    global ODE;    ODE = [];							% Custom integrator handle
    global STAGES; STAGES = 3;							% Default integrator configuration
    global RANK;   RANK = Inf;							% Maximum rank of decompositions

    persistent WC;
    persistent WO;
    persistent WQ;
    persistent WX;
    persistent FC;
    persistent FO;
    persistent FQ;
    persistent FX;
    persistent GC;
    persistent GO;
    persistent GQ;
    persistent GX;

    sysdim = [sys.M, sys.N, sys.Q];						% System dimension
    timdis = [sys.dt, sys.Tf];							% Time discretizations

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DEFAULT VALUES

    pr = hasfield(sys,'p',[]);                              			% Parameters
    nf = zeros(1,13);								% Configuration Flags
    ut = [];									% Training Input
    us = hasfield(sys,'us',zeros(sys.M,1));					% Steady-State Input
    xs = hasfield(sys,'xs',zeros(sys.N,1));					% Steady-State
    um = [];									% Input Perturbation Scales
    xm = [];									% Steady-State Perturbation Scales

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DECODE CONFIGURATION

    % Choose system integrator
    if isfield(config,'solver')

        if isa(config.solver,'function_handle')

            ODE = config.solver;						% Custom function handle
        else

            switch lower(config.solver)

                case 'rk1ex',   STAGES = 1;					% Explicit 1st Order Runge-Kutta Method (Explicit Euler's Method)

                case 'rk2ex',   STAGES = 2;					% Explicit 2nd Order Runge-Kutta Method (Heun's Method) Optimal Strong Stabilty Preserving, Low-Storage

                case 'rk45ex',  ODE = @rk45ex;					% Adaptive Embedded 4th/5th Order Runge-Kutta Method (Dormand-Prince Method)
            end%switch
        end%if
    end%if

    % Choose gramian kernel
    if isfield(config,'kernel') && isa(config.kernel,'function_handle')

        dp = config.kernel;
    end%if

    gtimes = @(m) m * m';

    dp = match(config,'kernel',[],{'sum',           @kernel_sum; ...		% Sum Peudo Kernel
                                   'trace',         @kernel_trace; ...		% Trace Peudo Kernel
                                   'diagonal',      @kernel_diagonal; ...	% Diagonal Pseudo Kernel
                                   'dmd',           @dmd; ...			% DMD Pseudo Kernel
                                   'position',      @(x,y) x(1:size(x,1)/2,:) * y(:,1:size(y,2)/2); ...		% Position Pseudo Kernel
                                   'velocity',      @(x,y) x(size(x,1)/2+1:end,:) * y(:,size(y,2)/2+1:end); ...	% Velocity Pseudo Kernel
                                   'quadratic',     @(x,y) (x * y).^2 + 1.0; ...		% Quadratic (Polynomial) Kernel
                                   'cubic',         @(x,y) (x * y).^3 + 1.0; ...		% Cubic (Polynomial) Kernel
                                   'sigmoid',       @(x,y) tanh(x * y - 1.0); ...		% Sigmoid Kernel
                                   'mercersigmoid', @(x,y) tanh(x - 1.0) * tanh(y - 1.0); ...	% Sigmoid-Mercer Kernel
                                   'logarithmic',   @(x,y) log(x + 1.0) * log(y + 1.0); ...	% Logarithmic Kernel
                                   'exponential',   @(x,y) exp(x * y); ...			% Exponential Kernel
                                   'gauss',         @(x,y) exp(-0.5 * gtimes(x-y')); ...	% Gauss Kernel
                                   'single',        @(x,y) single(x) * single(y)});		% Single Precision Kernel

    % Choose training input
    ut = match(config,'training','i',{'impulse', 'i'; ...			% Impulse input
                                      'step',    's'; ...			% Step input
                                      'chirp',   'h'; ...			% Chirp input
                                      'sinc',    'a'; ...			% Sinc input
                                      'random',  'r'});			% Random-binary input

    % Choose trajectory weighting
    nf(13) = match(config,'weighting',0,{'none',       0; ...                  % No weighting
                                         'linear',     1; ...			% Linear Time-Weighting
                                         'quadratic',  2; ...			% Quadratic Time-Weighting
                                         'state',      3; ...			% State-Based Weighting
                                         'scale',      4});			% Scale-Based Weighting

    % Choose trajectory centering
    nf(1) = match(config,'centering',0,{'none',     0; ...			% No Centering
                                        'steady',   1; ...			% Steady State
                                        'final',    2; ...			% Final State
                                        'mean',     3; ...			% Arithmetic Mean
                                        'rms',      4; ...			% Root-Mean-Squared
                                        'midrange', 5});			% Mid-Range

    % Choose perturbation scales
    nf([2,3]) = match(config,'scales',0,{'single',      0; ...			% No Subdivision: [1.0]
                                         'linear',      1; ...			% Linear Scale Subdivision: [0.25, 0.5, 0.75, 1.0]
                                         'geometric',   2; ...			% Geometric Scale Subdivision: [0.125, 0.25, 0.5, 1.0]
                                         'logarithmic', 3; ...			% Logarithmic Scale Subdivision: [0.001, 0.01, 0.1, 1.0]
                                         'sparse',      4});			% Sparse Scale Subdivision: [0.01, 0.5, 0.99, 1.0]

    % Choose perturbation rotations
    nf([4,5]) = match(config,'rotations',0,{'posneg', 0; ...			% Positive and negative rotations
                                            'single', 1});			% Only Positive Perturbations

    % Choose gramian normalization
    nf(6) = match(config,'normalize',0,{'none',   0; ...			% No normalization
                                        'steady', 1; ...			% Steady-State Normalization
                                        'jacobi', 2});				% Jacobi Normalization

    % State gramian variant
    nf(7) = match(config,'stype',0,{'standard',                0; ...		% Regular state Gramian
                                    'special' ,                1; ...		% Generic Non-Standard
                                    'output_controllability',  1; ...		% Output Controllabilty Gramian
                                    'averaged_observability',  1; ...		% Averaged Observability Gramian
                                    'nonsymmetric_minimality', 1});		% Nonsymmetric Cross Gramian

    % Extra input for observability and sensitivity
    nf(8) = match(config,'extra_input',0,{'none', 0; ...			% No extra input
                                          'yes',  1});				% Use extra input

    % Choose parameter centering
    nf(9) = match(config,'pcentering',0,{'none',        0; ...			% No Scaling
                                         'linear',      1; ...			% Linear Scaling and Parameter Centering
                                         'logarithmic', 2});			% Logarithmic Scaling Parameter Centering

    % Parameter gramian variant
    nf(10) = match(config,'ptype',0,{'standard',       0; ...			% Regular parameter Gramian
                                     'special',        1; ...			% Generic non-standard
                                     'io_sensitivity', 1; ...			% input-output-based sensitivity gramian
                                     'coarse_schur',   1});			% (cross-)identifiability gramian via approximate schur complement

    % Set maximum rank for decompositions
    RANK = hasfield(config,'max_order',Inf);

    % rom_training

    islinear = isfield(config,'linearity') && isequal(config.linearity,'linear');

    if islinear

        f = {sys.f, sys.F, sys.f};
        g = {sys.g, 1, sys.F};
        w = {'c', 'c', 'y'};
    else

        f = {sys.f, sys.f, sys.f};
        g = {sys.g, sys.g, sys.g};
        w = {'c', 'o', 'x'};
    end%if

    switch lower(task.type)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MATRIX EQUATIONS

        case 'matrix_equation' % OK

            v = match(task,'method',[],{'lyapunov',  1; ...
                                        'sylvester', 3});

            assert(not(isempty(v)),'est: Unknown matrix_equation method');

            r = emgr(f{v},g{v},sysdim,timdis,w{v},pr,nf,ut,us,xs,um,xm,dp);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SINGULAR VALUES

        case 'singular_values' % OK

            v = match(task,'method',[],{'controllability', 1; ...
                                        'observability',   2; ...
                                        'minimality',      3});

            assert(not(isempty(v)),'est: Unknown singular_value method');

            r = SVD(emgr(f{v},g{v},sysdim,timdis,w{v},pr,nf,ut,us,xs,um,xm,dp));

            if hasfield(config,'score',false)

                r = morscore(1:numel(r), r ./ max(r));
            elseif nargout == 2

                m = morscore(1:numel(r), r ./ max(r));
            end%if

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MODEL REDUCTION

        case 'model_reduction' % OK

            if isequal(task.method,'dmd_galerkin')

                dp = @dmd;
            end%if

            if     (isequal(task.method,'dmd_galerkin') || isequal(task.method,'poor_man')) && isequal(task.variant,'observability')

                W = {emgr(f{2},g{2},sysdim,timdis,w{2},pr,nf,ut,us,xs,um,xm,dp)};

            elseif (isequal(task.method,'dmd_galerkin') || isequal(task.method,'poor_man')) 

                W = {emgr(f{1},g{1},sysdim,timdis,w{1},pr,nf,ut,us,xs,um,xm,dp)};

            elseif isequal(task.variant,'observability')

                W = {emgr(f{1},g{1},sysdim,timdis,w{1},pr,nf,ut,us,xs,um,xm,dp), ...
                     emgr(f{2},g{2},sysdim,timdis,w{2},pr,nf,ut,us,xs,um,xm,dp)};

            elseif isequal(task.variant,'minimality')

                W = {emgr(f{3},g{3},sysdim,timdis,w{3},pr,nf,ut,us,xs,um,xm,dp)};
            else

                error('est: Unknown model_reduction variant');
            end%if

            reductor = match(task,'method',[],{'poor_man',            @poor_man;  ...
                                               'dmd_galerkin',        @poor_man;  ...
                                               'dominant_subspaces',  @dominant_subspaces; ...
                                               'approx_balancing',    @approx_balancing; ...
                                               'balanced_pod',        @balanced_pod; ...
                                               'balanced_truncation', @balanced_truncation}); 

            assert(not(isempty(reductor)),'est: Unknown model_reduction method');

            [UX,~,VX] = reductor(W);

            if any(strcmp(hasfield(config,'kernel',''),{'position','velocity'}))

                UX = blkdiag(UX,UX);
                VX = blkdiag(VX,VX);
            end%if

            if hasfield(config,'test',false)

                r = assess(sys,config,UX,VX,1,1);

                if hasfield(config,'score',false)

                    r = morscore(r{1},r{4});
                elseif nargout == 2

                    m = morscore(r{1},r{4});
                end%if
            else

                r = {UX,VX};
            end%if

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PARAMETER REDUCTION

        case 'parameter_reduction' % OK

            w = match(task,'method',[],{'observability', 'i'; ...
                                        'minimality',    'j'});

            assert(not(isempty(w)),'est: Unknown parameter_reduction method');

            W = emgr(sys.f,sys.g,sysdim,timdis,w,pr,nf,ut,us,xs,um,xm,dp);

            [UP,~,~] = SVD(W{2});

            if hasfield(config,'test',false)

                r = assess(sys,config,1,1,UP,UP);

                if hasfield(config,'score',false)

                    r = morscore(r{2},r{4});
                elseif nargout == 2

                    m = morscore(r{2},r{4});
                end%if
            else

                r = {UP,UP};
            end%if

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% COMBINED REDUCTION

        case 'combined_reduction' % OK

            if isequal(task.method,'observability')

                W = [emgr(sys.f,sys.g,sysdim,timdis,'c',pr,nf,ut,us,xs,um,xm,dp), ...
                     emgr(sys.f,sys.g,sysdim,timdis,'i',pr,nf,ut,us,xs,um,xm,dp)];

            elseif isequal(task.method,'minimality')

                W = emgr(sys.f,sys.g,sysdim,timdis,'j',pr,nf,ut,us,xs,um,xm,dp);
            else

                error('est: Unknown combined_reduction method');
            end%if

            [UP,~,~] = SVD(W{end});

            reductor = match(task,'variant',[],{'poor_man',            @poor_man;  ...
                                                'dominant_subspaces',  @dominant_subspaces; ...
                                                'approx_balancing',    @approx_balancing; ...
                                                'balanced_pod',        @balanced_pod; ...
                                                'balanced_truncation', @balanced_truncation});

            assert(not(isempty(reductor)),'est: Unknown combined_reduction variant');

            [UX,~,VX] = reductor(W(1:end-1));

            if hasfield(config,'test',false)

                r = assess(sys,config,UX,VX,UP,UP);

                if hasfield(config,'score',false)

                    r = morscore({r{1},r{2}},r{4});
                elseif nargout == 2

                    m = morscore({r{1},r{2}},r{4});
                end%if
            else

                r = {UX, VX; ...
                     UP, UP};
            end%if

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DECENTRALIZED CONTROL

        case 'decentralized_control' % OK

            nf(7) = 1;

            elem = @(v,k) v(k);

            coher = @(m) trace(m)^2 / sum(sum(m .* m'));

            gtrace = @(m) sum(sum(m .* m'));

            ys = sys.g(xs,us,pr,0);

            if islinear

                eg = @(ui,yj,dp) emgr(@(x,u,p,t) sys.f(x,us + sparse(ui,1,u,sys.M,1),p,t), ...
                                      @(x,u,p,t) sys.F(x,ys + sparse(yj,1,u,sys.Q,1),p,t), ...
                                      [1,sys.N,1],timdis,'y',pr,nf,ut,[],xs,um,xm,dp);
            else

                eg = @(ui,yj,dp) emgr(@(x,u,p,t) sys.f(x,us + sparse(ui,1,u,sys.M,1),p,t), ...
                                      @(x,u,p,t) elem(sys.g(x,u,p,t),yj), ...
                                      [1,sys.N,1],timdis,'x',pr,nf,ut,[],xs,um,xm,dp);
            end%if

            em = match(task,'method',[],{'relative_gain_array',  @(ui,yj) eg(ui,yj,@kernel_trace); ...
                                         'io_coherence',         @(ui,yj) coher(eg(ui,yj,[])); ...
                                         'io_pairing',           @(ui,yj) abs(det(eg(ui,yj,[]))); ...
                                         'participation_matrix', @(ui,yj) sqrt(gtrace(eg(ui,yj,[]))); ...
                                         'hardy_2',              @(ui,yj) abs(emgr(@(x,u,p,t) sys.f(x,us + sparse(ui,1,u,sys.M,1),p,t), ...
                                                                                   @(x,u,p,t) elem(sys.g(x,u,p,t),yj), ...
                                                                                   [1,sys.N,1],timdis,'c',pr,nf,ut,[],xs,um,xm)); ...
                                         'hardy_inf',            @(ui,yj) sum(abs(EIG(eg(ui,yj,[])))); ...
                                         'hankel_interaction',   @(ui,yj) abs(eigs(eg(ui,yj,[]),1)); ...
                                         'rms_hsv',              @(ui,yj) sum(SVD(eg(ui,yj,[])).^4)});

            assert(not(isempty(em)),'est: Unknown decentralized_control method');

            r = arrayfun(em,repmat(1:sys.Q,sys.M,1),repmat((1:sys.M)',1,sys.Q));
            if isequal(task.method,'relative_gain_array'), r = r.*pinv(r)'; end%if
            r = r./max(r(:));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% STATE SENSITIVITY

        case 'state_sensitivity' % OK

            v = match(task,'method',[],{'controllability', 1; ...
                                        'observability',   2; ...
                                        'minimality',      3});

            assert(not(isempty(v)),'est: Unknown state_sensitivity method');

            r = sqrt(abs(emgr(f{v},g{v},sysdim,timdis,w{v},pr,nf,ut,us,xs,um,xm,@kernel_diagonal)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PARAMETER SENSITIVITY

        case 'parameter_sensitivity' % OK

            v = match(task,'method',[],{'controllability', 0; ...
                                        'observability',   0; ...
                                        'minimality',      1});

            assert(not(isempty(v)),'est: Unknown parameter_sensitivity method');

            nf(10) = v;

            nf(7) = match(task,'method',0,{'observability', 1});

            ws = emgr(sys.f,sys.g,sysdim,timdis,'s',pr,nf,ut,us,xs,um,xm,dp);

            r = ws{2};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PARAMETER IDENTIFIABILITY

        case 'parameter_identifiability' % OK

            w = match(task,'method',[],{'observability', 'i'; ...
                                        'minimality',    'j'});

            assert(not(isempty(w)),'est: Unknown parameter_identifiability method');

            wi = emgr(sys.f,sys.g,sysdim,timdis,w,pr,nf,ut,us,xs,um,xm,dp);

            r = SVD(wi{2});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% UNCERTAINTY QUANTIFICATION

        case 'uncertainty_quantification' % OK

            v = match(task,'method',0,{'controllability', 0; ...
                                       'observability',   1});

            assert(not(isempty(w)),'est: Unknown uncertainty_quantification method');

            nf(7) = v;

            [UC,SC,~] = SVD(emgr(sys.f,sys.g,sysdim,timdis,'c',pr,nf,ut,us,xs,um,xm,dp));

            r = SC;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% NONLINEARITY QUANTIFICATION

        case 'nonlinearity_quantification' % OK

            w = match(task,'method',[],{'controllability', 'c'; ...
                                        'observability',   'o'; ...
                                        'minimality',      'x'; ...
                                        'correlation',     '!'});

            assert(not(isempty(w)),'est: Unknown nonlinearity_quantification method');

            if isequal(w,'!')

                rc = est(sys,setfield(task,'method','controllability'),config);
                ro = est(sys,setfield(task,'method','observability'),config);
                rx = est(sys,setfield(task,'method','minimality'),config);

                r = (rx .* rx) ./ (rc .* ro);
            else

                r = arrayfun(@(k) emgr(sys.f,sys.g,sysdim,timdis,w,pr,nf,ut,us,xs,k,k,@kernel_trace),linspace(1.0,10.0,10));
            end%if

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GRAMIAN INDEX

        case 'gramian_index' % OK

            inw = match(task,'method',[],{'sigma_min',             @(w) min(SVD(w)); ...
...
                                          'harmonic_mean',         @(w) size(w,1)/sum(1./SVD(w))
...
                                          'geometric_mean',        @(w) prod(SVD(w))^(1.0/size(w,1)); ...
...
                                          'energy_fraction',       @(w) sum(SVD(w)); ...
...
                                          'operator_norm',         @(w) norm(w,'fro'); ...
...
                                          'sigma_max',             @(w) svds(w,1); ...
...
                                          'log_det',               @(w) log(prod(SVD(w))); ...
...
                                          'storage_efficiency',    @(w) sqrt(prod(SVD(w))/prod(diag(w))); ...
...
                                          'unobservability_index', @(w) 1.0./sqrt(min(SVD(w))); ...
...
                                          'performance_index',     @(w) trace(w) * prod(SVD(w))^(1.0/size(w,1))});

            assert(not(isempty(inw)),'est: Unknown gramian_index method');

            switch hasfield(task,'variant','')

                case 'controllability'

                    if isempty(WC) || not(isequal(FC,sys.f)) || not(isequal(GC,sys.g))

                        WC = emgr(f{1},g{1},sysdim,timdis,w{1},pr,nf,ut,us,xs,um,xm,dp);
                        FC = sys.f;
                        GC = sys.g;
                    end%if

                    W = WC;

                case 'observability'

                    if isempty(WO) || not(isequal(FO,sys.f)) || not(isequal(GO,sys.g))

                        WO = emgr(f{2},g{2},sysdim,timdis,w{2},pr,nf,ut,us,xs,um,xm,dp);
                        FO = sys.f;
                        GO = sys.g;
                    end%if

                    W = WC;

                case 'minimality'

                    if isempty(WX) || not(isequal(FX,sys.f)) || not(isequal(GX,sys.g))

                        WX = emgr(f{3},g{3},sysdim,timdis,w{3},pr,nf,ut,us,xs,um,xm,dp);
                        FX = sys.f;
                        GX = sys.g;
                    end%if

                    W = WX;
            end%switch

            r = eps + abs(inw(W) - arrayfun(@(k) inw(sys.proj{2}(:,1:k)'*W*sys.proj{1}(:,1:k)),1:min(sys.N,RANK)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SYSTEM INDEX

        case 'system_index' % OK

            inwx = match(task,'method',[],{'cauchy_index',    @(wx) sum(sign(real(EIG(wx)))); ...
...
                                           'system_entropy',  @(wx) size(wx,1)/log(2.0*exp(1.0)*pi) + sum(log(abs(abs(EIG(wx))))); ...
...
                                           'system_symmetry', @(wx) sqrt(abs(sum(sum(wx.*wx'))))/norm(wx,'Fro'); ...
...
                                           'io_coherence',    @(wx) abs(sum(sum(wx.*wx')))/trace(wx)^2; ...
...
                                           'system_gain',     @(wx) abs(trace(wx))});

            if not(isempty(inwx))

                if isempty(WX) || not(isequal(FX,sys.f)) || not(isequal(GX,sys.g))

                    WX = emgr(f{3},g{3},sysdim,timdis,w{3},pr,nf,ut,us,xs,um,xm,dp);
                    FX = sys.f;
                    GX = sys.g;
                end%if

                r = eps + abs(inwx(WX) - arrayfun(@(k) inwx(sys.proj{2}(:,1:k)'*WX*sys.proj{1}(:,1:k)),1:min(sys.N,RANK)));
                return
            end%if

            inco = match(task,'method',[],{'gramian_distance',            @(wc,wo) norm(log(sqrt(EIG(wc*wo))),2); ...
...
                                           'network_sensitivity',         @(wc,wo) trace(wc) + trace(wo); ...
...
                                           'geometric_mean_hsv',          @(wc,wo) prod(sqrt(EIG(wc*wo)))^(1.0/size(wc,1)); ...
...
                                           'rv_coefficient',              @(wc,wo) sum(sum(wc.*wo))/(norm(wc,'Fro')*norm(wo,'Fro'))});

            if not(isempty(inco))

                if isempty(WC) || not(isequal(FC,sys.f)) || not(isequal(GC,sys.g)) || not(isequal(FO,sys.f)) || not(isequal(GO,sys.g))

                    nf(7) = 0;
                    WC = emgr(f{1},g{1},sysdim,timdis,w{1},pr,nf,ut,us,xs,um,xm,dp);
                    WO = emgr(f{2},g{2},sysdim,timdis,w{2},pr,nf,ut,us,xs,um,xm,dp);
                    FC = sys.f;
                    GC = sys.g;
                    FO = sys.f;
                    GO = sys.g;
                end%if

                r = eps + abs(inco(WC,WO) - arrayfun(@(k) inco(sys.proj{2}(:,1:k)'*WC*sys.proj{1}(:,1:k),sys.proj{2}(:,1:k)'*WO*sys.proj{1}(:,1:k)),1:min(sys.N,RANK)));
                return
            end%if

            error('est: Unknown system_index method');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SYSTEM NORM

        case 'system_norm' % OK

            inoc = match(task,'method',[],{'hardy_2_norm', @(w) sqrt(abs(trace(w)))});

            if not(isempty(inoc))

                if isempty(WQ) || not(isequal(FQ,sys.f)) || not(isequal(GQ,sys.g))

                    nf(7) = 1;
                    WQ = emgr(f{1},g{1},sysdim,timdis,w{1},pr,nf,ut,us,xs,um,xm,dp);
                    FQ = sys.f;
                    GQ = sys.g;
                end%if

                r = eps + abs(inoc(WQ) - arrayfun(@(k) inoc(emgr(f{1},@(x,u,p,t) g{1}(sys.proj{1}(:,1:k)*(sys.proj{2}(:,1:k)'*x),u,p,t),sysdim,timdis,w{1},pr,nf,ut,us,xs,um,xm,dp)),1:min(sys.N,RANK))); 
                return
            end%if

            inco = match(task,'method',[],{'hardy_inf_norm',              @(wc,wo) sum(sqrt(EIG(wc*wo))); ...
...
                                           'hilbert_schmidt_hankel_norm', @(wc,wo) norm(wc*wo,'Fro'); ...
...
                                           'hankel_norm',                 @(wc,wo) sqrt(min(EIG(wc*wo)))});

            if not(isempty(inco))

                if isempty(WC) || not(isequal(FC,sys.f)) || not(isequal(GC,sys.g)) || not(isequal(FO,sys.f)) || not(isequal(GO,sys.g))

                    nf(7) = 0;
                    WC = emgr(f{1},g{1},sysdim,timdis,w{1},pr,nf,ut,us,xs,um,xm,dp);
                    WO = emgr(f{2},g{2},sysdim,timdis,w{2},pr,nf,ut,us,xs,um,xm,dp);
                    FC = sys.f;
                    GC = sys.g;
                    FO = sys.f;
                    GO = sys.g;
                end%if

                r = eps + abs(inco(WC,WO) - arrayfun(@(k) inco(sys.proj{2}(:,1:k)'*WC*sys.proj{1}(:,1:k),sys.proj{2}(:,1:k)'*WO*sys.proj{1}(:,1:k)),1:min(sys.N,RANK)));
                return
            end%if

            error('est: Unknown system_norm method');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TAU FUNCTION

        case 'tau_function'

            r = arrayfun(@(k) prod(real(EIG(eye(sys.N) + emgr(f{3},g{3},sysdim,timdis,w{3},pr,nf,ut,us,xs,um,xm,@(x,y) x(:,k:end)*y(k:end,:))))),1:floor(sys.Tf / sys.dt));
    end%switch

    ODE = [];
    STAGES = [];
    RANK = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DECOMPOSITION WRAPPER

function varargout = SVD(A)
% summary: svd/svds wrapper

    global RANK;

    switch nargout

        case 1

            if isinf(RANK)
                varargout = {svd(A)};
            else
                varargout = {svds(A,r)};
            end%if

        case 2

            if isinf(RANK)
                [U,D,~] = svd(A);
            else
                [U,D,~] = svds(A,r);
            end%if
            D = diag(D);
            varargout = {U,D};

        case 3

            if isinf(RANK)
                [U,D,V] = svd(A);
            else
                [U,D,V] = svds(A,r);
            end%if
            D = diag(D);
            varargout = {U,D,V};

    end%switch
end

function varargout = EIG(A)
% summary: eig/eigs wrapper

    global RANK;

    switch nargout

        case 1

            if isinf(RANK)
                varargout = {eig(A)};
            else
                varargout = {eigs(A,r)};
            end%if

        case 2

            if isinf(RANK)
                [U,D,~] = eig(A,'vector');
            else
                [U,D,~] = eigs(A,r);
                D = diag(D);
            end%if
            varargout = {U,D};

    end%switch
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% UTILITIES

function r = hasfield(str,key,def)
% summary: get field key from struct str if exists otherwise return def

    if isfield(str,key)

        r = getfield(str,key);
    else

        r = def;
    end%if
end

function r = match(str,key,def,map)
% summary: return map(ped) values for member key of struct str otherwise def

    s = cell2struct(map(:,2),map(:,1),1);

    if not(isfield(str,key)) || not(isfield(s,getfield(str,key)))

        r = def;
    else

        r = getfield(s,getfield(str,key));
    end%if
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PSEUDO-KERNELS

function r = kernel_sum(x,y)
% summary: Sum Pseudo-Kernel

    r = sum(sum(x*y));
end

function r = kernel_trace(x,y)
% summary: Trace Pseudo-Kernel

    r = sum(sum(x.*y'));
end

function r = kernel_diagonal(x,y)
% summary: Diagonal Pseudo-Kernel

    r = sum(x.*y',2);
end

function r = dmd(x,y)
% summary: Dynamic-Mode-Decomposition-Galerkin Pseudo Kernel

    r = x(:,2:end) * pinv(y(1:end-1,:)',sqrt(eps));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ADAPTIVE RUNGE-KUTTA SOLVER

function y = rk45ex(f,g,t,x0,u,p)
% summary: Adaptive 4th/5th Dormand-Prince Runge-Kutta method

    [S,x] = ode45(@(t,x) f(x,u(t),p,t),[0,t(2)],x0,'InitialStep',t(1));
    z = arrayfun(@(k) g(x(k,:)',u(S(k)),p,S(k)),1:numel(S),'UniformOutput',false);
    y = interp1(S,z',0:t(1):t(2))';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MORSCORE (MODEL ORDER REDUCTION SCORE)

function s = morscore(orders,errors)

    if iscell(orders) && all(size(errors)>1)

        nx = orders{1} ./ max(orders{1});
        ny = orders{2} ./ max(orders{2});
        nz = log10(errors + eps) ./ floor(log10(eps));

        s = max(0,trapz(ny(:),trapz(nx(:),nz,2)));
    else

        nx = orders ./ max(orders);
        ny = log10(errors + eps) ./ floor(log10(eps));

        s = max(0,trapz(nx(:),ny(:)));
    end%if
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% REDUCTORS

function [U,D,V] = poor_man(W)
% summary: Poor man's method (pod)

    [U,D,~] = SVD(W{1});
    V = U;
end

function [U,D,V] = dominant_subspaces(W)
% summary: Dominant subspaces

    if isequal(numel(W),1)

        [UX,DX,VX] = SVD(W{1});
        [U,D,~] = SVD([UX.*DX',VX.*DX']);
    else

        [UC,DC,~] = SVD(W{1});
        [UO,DO,~] = SVD(W{2});
        [U,D,~] = SVD([UC.*DC',UO.*DO']);
    end%if

    V = U;
end

function [U,D,V] = approx_balancing(W)
% summary: approximate balancing (modified pod)

    if isequal(numel(W),1)

        [U,DX,VX] = SVD(W{1});
        D = diag(DX);
        V = VX*(VX'*U);
    else

        [U,DC,~] = SVD(W{1});
        [VX,DO,~] = SVD(W{2});
        D = diag(DC)./diag(DO);
        V = VX*(VX'*U);
    end%if
end

function [U,D,V] = balanced_pod(W)
% summary: Balanced pod

    if isequal(numel(W),1)

        [LC,EC] = SVD(W{1});
        [LO,EO] = SVD(W{1}');
    else

        [LC,EC] = SVD(W{1});
        [LO,EO] = SVD(W{2});
    end%if

    LC = LC .* sqrt(abs(EC))';
    LO = LO .* sqrt(abs(EO))';

    [UB,HSV,VB] = svd(LC' * LO,'econ');
    D = sqrt(diag(HSV) + 2.0*eps)';
    U = LC * (UB ./ D);
    V = LO * (VB ./ D);
end

function [U,D,V] = balanced_truncation(W)
% summary: Balanced truncation

    if isequal(numel(W),1)

        [LC,EC] = EIG(W{1});
        [LO,EO] = EIG(W{1}');
    else

        [LC,EC] = EIG(W{1}*W{2});
        [LO,EO] = EIG(W{2}*W{1});
    end%if

    LC = LC .* sqrt(abs(EC))';
    LO = LO .* sqrt(abs(EO))';

    [UB,HSV,VB] = svd(LC' * LO,'econ');
    D = sqrt(diag(HSV) + 2.0*eps)';
    U = LC * (UB ./ D);
    V = LO * (VB ./ D);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% REDUCED ORDER MODEL EVALUATION

function r = assess(sys,config,XL,XR,PL,PR)
% summary:

    global ODE;

    pr = hasfield(sys,'p',0);
    us = hasfield(sys,'us',zeros(sys.M,1));
    xs = hasfield(sys,'xs',zeros(sys.N,1));
    x0 = hasfield(sys,'x0',zeros(sys.N,1));

    rand('seed',1009);
    ur = rand(1,floor(sys.Tf / sys.dt) + 1);
    u = @(t) ur(1,floor(t / sys.dt) + 1);

    skip_x = hasfield(config,'skip_x',1);
    skip_p = hasfield(config,'skip_p',1);
    num_test_param = hasfield(config,'num_test_param',1);

    if isequal(num_test_param,1) || isequal(size(pr,2),1)
        param = pr;
    else
        pmin = min(pr,[],2);
        pmax = max(pr,[],2);
        param = pmin + abs(pmax - pmin) .* rand(size(pr,1),num_test_param);
    end%if

    if isequal(numel(XL),1)
        max_x = 1;
    else
        max_x = min(size(XL,1)-1,size(XL,2));
    end%if

    if isequal(numel(PL),1)
        max_p = 1;
    else
        max_p = min(size(PL,1)-1,size(PL,2));
    end%if

    test_x = skip_x:skip_x:max_x;
    test_p = skip_p:skip_p:max_p;

    norms = { @(y) sys.dt * norm(y(:),1), ...					% L1 time series norm
              @(y) sqrt(sys.dt) * norm(y(:),2), ...				% L2 time series norm
              @(y) norm(y(:),Inf), ...						% Linf time series norm
              @(y) sum(abs(prod(y,1).^(1/size(y,1))))};				% L0 time series norm

    ln = cellfun(@(n) zeros(numel(test_x),numel(test_p)),norms,'UniformOutput',false);

    for q = 1:num_test_param

        Y = ODE(sys.f,sys.g,[sys.dt,sys.Tf],x0,@(t) us + u(t),param(:,q));

        for n = test_x

            xl = XL(:,1:n);
            xr = XR(:,1:n)';

            ix = find(test_x==n);

            for p = test_p

                pl = PL(:,1:p);
                pr = PR(:,1:p)';

                ip = find(test_p==p);

                y = ODE(@(x,u,p,t) xr*sys.f(xs + xl*x,u,p,t), ...
                        @(x,u,p,t)    sys.g(xs + xl*x,u,p,t), ...
                        [sys.dt,sys.Tf], xr*x0, @(t) us + u(t), pl*(pr*param(:,q)));

                for m = 1:numel(norms)

                    ln{m}(ix,ip) = ln{m}(ix,ip) + (norms{m}(Y - y)).^2;
                end%for
            end%for
        end%for
    end%for

    ln = cellfun(@(l) sqrt(l)./sqrt(max(l(:))),ln,'UniformOutput',false);

    r = [{test_x, test_p}, ln];
end
