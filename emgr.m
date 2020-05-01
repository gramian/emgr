function W = emgr(f,g,s,t,w,pr,nf,ut,us,xs,um,xm,dp)
%% emgr - EMpirical GRamian Framework
%
%  project: emgr ( https://gramian.de )
%  version: 5.8 (2020-05-01)
%  authors: C. Himpe (0000-0003-2194-6754)
%  license: BSD-2-Clause (opensource.org/licenses/BSD-2-Clause)
%  summary: Empirical system Gramians for (nonlinear) input-output systems.
%
% DESCRIPTION:
%
%  Empirical gramian matrix and empirical covariance matrix computation
%  for model reduction, decentralized control, nonlinearity quantification,
%  sensitivity analysis, parameter identification, uncertainty quantification &
%  combined state and parameter reduction of large-scale input-output systems.
%  Data-driven analysis of input-output coherence and system-gramian-based
%  nonlinear model order reduction. Compatible with OCTAVE and MATLAB.
%
% ALGORITHM:
%
%  C. Himpe (2018). emgr - The Empirical Gramian Framework. Algorithms 11(7):91
%  <https://doi.org/10.3390/a11070091 doi:10.3390/a11070091>
%
% USAGE:
%
%  W = emgr(f,g,s,t,w,[pr],[nf],[ut],[us],[xs],[um],[xm],[dp])
%
% MANDATORY ARGUMENTS:
%
%   f {handle} vector field: x' = f(x,u,p,t)
%   g {handle} output functional: y = g(x,u,p,t)
%   s {vector} system dimensions: [inputs, states, outputs]
%   t {vector} time discretization: [time-step, time-horizon]
%   w  {char}  empirical gramian type:
%    * 'c' empirical controllability gramian (Wc)
%    * 'o' empirical observability gramian (Wo)
%    * 'x' empirical cross gramian (Wx aka Wco)
%    * 'y' empirical linear cross gramian (Wy)
%    * 's' empirical sensitivity gramian (Ws)
%    * 'i' empirical identifiability gramian (Wi)
%    * 'j' empirical joint gramian (Wj)
%
% OPTIONAL ARGUMENTS:
%
%  pr {matrix|0} parameter vector(s), each column is one parameter sample
%  nf {vector|0} option flags, thirteen component vector, default all zero:
%    * centering: none(0), steady(1), last(2), mean(3), rms(4), midrange(5)
%    * input scales: single(0), linear(1), geometric(2), log(3), sparse(4)
%    * state scales: single(0), linear(1), geometric(2), log(3), sparse(4)
%    * input rotations: unit(0), single(1)
%    * state rotations: unit(0), single(1)
%    * normalization (only: Wc, Wo, Wx, Wy): none(0), steady(1), Jacobi(2)
%    * state gramian variant:
%      * controllability gramian type (only: Wc): regular(0), output(1)
%      * observability gramian type (only: Wo, Wi): regular(0), averaged(1)
%      * cross gramian type (only: Wx, Wy, Wj): regular(0), non-symmetric(1)
%    * extra input (only: Wo, Wx, Ws, Wi, Wj): none(0), yes(1)
%    * parameter centering (only: Ws, Wi, Wj): none(0), linear(1), log(2)
%    * parameter gramian variant:
%      * averaging type (only: Ws): input-state(0), input-output(1)
%      * Schur-complement (only: Wi, Wj): approx(0), coarse(1)
%    * cross gramian partition size (only: Wx, Wj): full(0), partitioned(<N)
%    * cross gramian partition index (only: Wx, Wj): partition(>0)
%    * weighting: none(0), time-linear(1), time-squared(2), state(3), scale(4)
%  ut {handle|'i'} input function: u_t = ut(t) or character:
%    * 'i' delta impulse input
%    * 's' step input / load vector / source term
%    * 'c' decaying exponential chirp input
%    * 'a' sinc (cardinal sine) input
%    * 'r' pseudo-random binary input
%  us {vector|0} steady-state input (1 or M rows)
%  xs {vector|0} steady-state and nominal initial state x_0 (1 or N rows)
%  um {matrix|1} input scales (1 or M rows)
%  xm {matrix|1} initial-state scales (1 or N rows)
%  dp {handle|@mtimes} inner product or kernel: xy = dp(x,y)
%
% RETURNS:
%
%  W {matrix} State-space Gramian Matrix (for: Wc, Wo, Wx, Wy)
%  W  {cell}  [State, Parameter]-space Gramian (for: Ws, Wi, Wj)
%
% CITE AS:
%
%  C. Himpe (2020). emgr - EMpirical GRamian Framework (Version 5.8)
%  [Software]. Available from https://gramian.de . doi:10.5281/zenodo.3779889
%
% KEYWORDS:
%
%  model reduction, system gramians, empirical gramians, cross gramian, MOR
%
% SEE ALSO: gram (Control System Toolbox)
%
% COPYRIGHT: Christian Himpe
%
% For more information, see: <https://gramian.de>

    % Integrator Handle (for custom solvers)
    global ODE;
    if not(isa(ODE,'function_handle')), ODE = @ssp2; end%if

    % Version Info (use also to export integrator handle)
    if isequal(f,'version'), W = 5.8; return; end%if

    % Default Arguments
    if (nargin <  6) || isempty(pr), pr = 0.0; end%if
    if (nargin <  7) || isempty(nf), nf = 0.0; end%if
    if (nargin <  8) || isempty(ut), ut = 'i'; end%if
    if (nargin <  9) || isempty(us), us = 0.0; end%if
    if (nargin < 10) || isempty(xs), xs = 0.0; end%if
    if (nargin < 11) || isempty(um), um = 1.0; end%if
    if (nargin < 12) || isempty(xm), xm = 1.0; end%if
    if (nargin < 13) || isempty(dp), dp = @mtimes; end%if

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SETUP

    % System Dimensions
    M = s(1);                 % Number of inputs
    N = s(2);                 % Number of states
    Q = s(3);                 % Number of outputs
    P = size(pr,1);           % Dimension of parameter
    K = size(pr,2);           % Number of parameter-samples

    % Time Discretization
    dt = t(1);                % Time-step width
    Tf = t(2);                % Time horizon
    nt = floor(Tf / dt) + 1;  % Number of time-steps plus initial value

    % Force lower-case Gramian type
    w = lower(w);

    % Lazy Output Functional
    if isnumeric(g) && isequal(g,1), g = @id; Q = N; end%if

    % Pad Configuration Flag Vector
    nf(end + 1:13) = 0;

    % Built-in Input Functions
    if not(isa(ut,'function_handle'))

        switch lower(ut)
            case 's'   % Step Input
                ut = @(t) 1.0;

            case 'c'   % Chirp Input
                a0 = (2.0 * pi) / (4.0 * dt) * Tf / log(4.0 * (dt / Tf));
                b0 = (4.0 * (dt / Tf)) ^ (1.0 / Tf);
                ut = @(t) 0.5 * cos(a0 * (b0 ^ t - 1.0)) + 0.5;

            case 'a'   % Sinc Input
                ut = @(t) sin(t / dt) / ((t / dt) + (t==0));

            case 'r'   % Binary Input
                rt = randi([0,1],1,nt);
                ut = @(t) rt(:,floor(t / dt) + 1);

            otherwise  % Impulse Input
                ut = @(t) (t <= dt) / dt;
        end%switch
    end%if

    % Lazy Optional Arguments
    if isscalar(us), us = repmat(us,M,1); end%if
    if isscalar(xs), xs = repmat(xs,N,1); end%if
    if isscalar(um), um = repmat(um,M,1); end%if
    if isscalar(xm), xm = repmat(xm,N,1); end%if

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CONFIGURATION

    % Trajectory Weighting
    switch nf(13)
        case 1     % Linear time-weighting
            wei = @(m) sqrt(0:dt:Tf);

        case 2     % Quadratic time-weighting
            wei = @(m) (0:dt:Tf) * sqrt(0.5);

        case 3     % State-weighting
            wei = @(m) vecnorm(m,2,1);

        case 4     % Scale-weighting
            wei = @(m) 1.0 ./ max(eps,vecnorm(m,Inf,2));

        otherwise  % None
            wei = @(m) 1.0;
    end%switch

    % Trajectory Centering
    switch nf(1)
        case 1     % Steady state / output
            avg = @(m,s) s;

        case 2     % Final state / output
            avg = @(m,s) m(:,end);

        case 3     % Temporal mean of state / output
            avg = @(m,s) mean(m,2);

        case 4     % Temporal root-mean-square of state / output
            avg = @(m,s) sqrt(mean(m .* m,2));

        case 5     % Temporal mid-range of state / output
            avg = @(m,s) 0.5 * (max(m,[],2) + min(m,[],2));

        otherwise  % None
            avg = @(m,s) 0;
    end%switch

    % Gramian Normalization
    if nf(6)

        TX = xs;             % Steady-state preconditioner

        if isequal(nf(6),2)  % Jacobi-type preconditioner

            NF = nf;
            NF(6) = 0;
            if isequal(w,'c') || isequal(w,'s'), NF(7) = 0; end%if
            WN = w;
            if isequal(w,'s'), WN = 'c'; end%if
            if isequal(w,'i'), WN = 'o'; end%if
            if isequal(w,'j'), WN = 'x'; end%if
            PR = mean(pr,2);
            DP = @(x,y) sum(x .* y',2);  % Diagonal-only pseudo-kernel
            TX = sqrt(abs(emgr(f,g,s,t,WN,PR,NF,ut,us,xs,um,xm,DP)));
        end%if

        TX(abs(TX) < (8*eps)) = 1.0;
        if isequal(w,'y'), tx = TX; else, tx = 1.0; end%if
        f = @(x,u,p,t) f(TX .* x,u,p,t) ./ TX;
        g = @(x,u,p,t) g(TX .* x,u,p,t) ./ tx;
        xs = xs ./ TX;
        nf(6) = 0;
    end%if

    % Non-symmetric cross Gramian and average observability Gramian
    if nf(7), R = 1; else, R = Q; end%if

    % Extra input for control explicit observability
    if nf(8), up = @(t) us + ut(t); else, up = @(t) us; end%if

    % Scale Sampling
    if isequal(size(um,2),1), um = um * scales(nf(2),nf(4)); end%if
    if isequal(size(xm,2),1), vm = xm(1:Q) * scales(nf(2),nf(4)); end%if
    if isequal(size(xm,2),1), xm = xm * scales(nf(3),nf(5)); end%if

    A = size(xm,1);  % Number of total states (regular plus augmented)
    C = size(um,2);  % Number of input scales sets
    D = size(xm,2);  % Number of state scales sets

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GRAMIAN COMPUTATION

    switch w  % Empirical system gramian types

        % Common Layout:
        %   For each {parameter, scale, input/state/parameter component}:
        %     Perturb, simulate, weight, center, normalize, accumulate
        %   Parameter gramians call state gramians

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Empirical Controllability Gramian

        case 'c'

            W = 0;  % Initialize gramian variable
            for k = 1:K
                for c = 1:C
                    for m = find(um(:,c))'  % parfor
                        em = sparse(m,1,um(m,c),M + P,1);
                        umc = @(t) up(t) + ut(t) .* em(1:M);
                        pmc = pr(:,k) + em(M + 1:end);
                        if nf(7)  % Output controllability gramian
                            x = ODE(f,g,t,xs,umc,pmc);
                        else      % Regular controllability gramian
                            x = ODE(f,@id,t,xs,umc,pmc);
                        end%if
                        x = x .* wei(x);
                        x = x - avg(x,xs);
                        x = x ./ um(m,c);
                        W = W + dp(x,x');
                    end%for
                end%for
            end%for
            W = W * (dt / (C * K));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Empirical Observability Gramian

        case 'o'

            W = 0;               % Initialize gramian variable
            o = zeros(R*nt, A);  % Pre-allocate observability matrix
            for k = 1:K
                for d = 1:D
                    for n = find(xm(:,d))'  % parfor
                        en = sparse(n,1,xm(n,d),N + P,1);
                        xnd = xs + en(1:N);
                        pnd = pr(:,k) + en(N + 1:end);
                        y = ODE(f,g,t,xnd,up,pnd);
                        y = y .* wei(y);
                        y = y - avg(y,g(xs,us,pnd,0));
                        y = y ./ xm(n,d);
                        if nf(7)  % Average observability gramian
                            o(:,n) = sum(y,1)';
                        else      % Regular observability gramian
                            o(:,n) = y(:);
                        end%if
                    end%for
                    W = W + dp(o',o);
                end%for
            end%for
            W = W * (dt / (D * K));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Empirical Cross Gramian

        case 'x'

            assert(isequal(M,Q) || nf(7),' emgr: non-square system!');

            i0 = 1;  % Start column index
            i1 = A;  % Final column index

            % Partitioned cross gramian
            if nf(11) > 0
                sp = round(nf(11));  % Partition size
                ip = round(nf(12));  % Partition index
                i0 = i0 + (ip - 1) * sp;
                i1 = min(i0 + sp - 1,N);
                if i0 > N
                    i0 = i0 - (ceil(N / sp) * sp - N);
                    i1 = min(i0 + sp - 1,A);
                end%if

                if (ip < 1) || (i0 > i1) || (i0 < 0), W = 0; return; end%if
            end%if

            W = 0;                      % Initialize gramian variable
            o = zeros(nt, i1-i0+1, R);  % Pre-allocate observability 3-tensor
            for k = 1:K
                for d = 1:D
                    for n = find(xm(i0:i1,d))'
                        en = sparse(i0-1+n,1,xm(i0-1+n,d),N+P,1);
                        xnd = xs + en(1:N);
                        pnd = pr(:,k) + en(N + 1:end);
                        y = ODE(f,g,t,xnd,up,pnd);
                        y = y .* wei(y);
                        y = y - avg(y,g(xs,us,pnd,0));
                        y = y ./ xm(i0 - 1 + n,d);
                        if nf(7)  % Non-symmetric cross gramian
                            o(:,n,1) = sum(y,1)';
                        else      % Regular cross gramian
                            o(:,n,:) = y';
                        end%if
                    end%for
                    for c = 1:C  % parfor
                        for m = find(um(:,c))'
                            em = sparse(m,1,um(m,c),M,1);
                            umc = @(t) us + ut(t) .* em;
                            x = ODE(f,@id,t,xs,umc,pr(:,k));
                            x = x .* wei(x);
                            x = x - avg(x,xs);
                            x = x ./ um(m,c);
                            if nf(7)  % Non-symmetric cross gramian
                                W = W + dp(x,o(:,:,1));
                            else      % Regular cross gramian
                                W = W + dp(x,o(:,:,m));
                            end%if
                        end%for
                    end%for
                end%for
            end%for
            W = W * (dt / (C * D * K));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Empirical Linear Cross Gramian

        case 'y'

            assert(isequal(M,Q) || nf(7),' emgr: non-square system!');
            assert(isequal(C,size(vm,2)),' emgr: scale count mismatch!');

            W = 0;               % Initialize gramian variable
            a = cell(1,R);       % Initialize adjoint cache
            a{1} = zeros(N,nt);  % Pre-allocate accumulator
            for k = 1:K
                for c = 1:C
                    for q = find(vm(:,c))'
                        em = sparse(q,1,vm(q,c),Q,1);
                        vqc = @(t) us + ut(t) .* em;
                        z = ODE(g,@id,t,xs,vqc,pr(:,k));
                        z = z .* wei(z);
                        z = z - avg(z,xs);
                        z = z ./ vm(q,c);
                        if nf(7)  % Non-symmetric cross gramian
                            a{1} = a{1} + z;
                        else      % Regular cross gramian
                            a{q} = z;
                        end%if
                    end%for
                    for m = find(um(:,c))'  % parfor
                        em = sparse(m,1,um(m,c),M,1);
                        umc = @(t) us + ut(t) .* em;
                        x = ODE(f,@id,t,xs,umc,pr(:,k));
                        x = x .* wei(x);
                        x = x - avg(x,xs);
                        x = x ./ um(m,c);
                        if nf(7)  % Non-symmetric cross gramian
                            W = W + dp(x,a{1}');
                        else      % Regular cross gramian
                            W = W + dp(x,a{m}');
                        end%if
                    end%for
                end%for
            end%for
            W = W * (dt / (C * K));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Empirical Sensitivity Gramian

        case 's'

            % Controllability gramian
            [pr,pm] = pscales(pr,nf(9),C);
            W{1} = emgr(f,g,s,t,'c',pr,nf,ut,us,xs,um,xm,dp);

            if not(nf(10))  % Input-state sensitivty gramian
                DP = @(x,y) sum(sum(x .* y'));       % Trace pseudo-kernel
            else            % Input-output sensitivity gramian
                DP = @(x,y) sum(reshape(y,R,[]));    % Custom pseudo-kernel
                Y = emgr(f,g,s,t,'o',pr,nf,ut,us,xs,um,xm,DP);
                DP = @(x,y) abs(sum(y(:) .* Y(:)));  % Custom pseudo-kernel
            end%if

            V = zeros(P,1);
            for p = 1:P  % parfor
                pmp = sparse(M + p,1:C,pm(p,:),M + P,C);
                V(p) = emgr(f,g,s,t,'c',pr,nf,ut,us,xs,pmp,xm,DP);
            end%for
            W{2} = V; % (Diagonal) Sensitivity gramian

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Empirical Identifiability Gramian

        case 'i'

            % Augmented observability gramian
            [pr,pm] = pscales(pr,nf(9),D);
            V = emgr(f,g,s,t,'o',pr,nf,ut,us,xs,um,[xm;pm],dp);

            W{1} = V(1:N, 1:N);    % Observability gramian
            WM = V(1:N, N+1:N+P);  % Mixed block

            % Identifiability gramian
            if not(nf(10))  % Schur-complement via approximate inverse
                W{2} = V(N+1:N+P, N+1:N+P) - (WM' * ainv(W{1}) * WM);   
            else            % Coarse Schur-complement via zero
                W{2} = V(N+1:N+P, N+1:N+P);
            end%if

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Empirical Joint Gramian

        case 'j'

            % Joint gramian
            [pr,pm] = pscales(pr,nf(9),D);
            V = emgr(f,g,s,t,'x',pr,nf,ut,us,xs,um,[xm;pm],dp);

            % Joint gramian partition
            if nf(11), W = V; return; end%if            

            W{1} = V(1:N, 1:N);    % Cross gramian
            WM = V(1:N, N+1:N+P);  % Mixed block

            % Cross-identifiability gramian
            if not(nf(10))  % Schur-complement via approximate inverse
                W{2} = 0.5 * (WM' * ainv(W{1} + W{1}') * WM);
            else            % Coarse Schur-complement via identity
                W{2} = 0.5 * (WM' * WM);                
            end%if

        otherwise

            error(' emgr: unknown empirical gramian type!');
    end%switch
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LOCAL FUNCTION: scales

function s = scales(nf1,nf2)
% summary: Input and initial state perturbation scales

    switch nf1
        case 1     % Linear
            s = [0.25, 0.50, 0.75, 1.0];

        case 2     % Geometric
            s = [0.125, 0.25, 0.5, 1.0];

        case 3     % Logarithmic
            s = [0.001, 0.01, 0.1, 1.0];

        case 4     % Sparse
            s = [0.01, 0.50, 0.99, 1.0];

        otherwise  % One
            s = 1.0;
    end%switch

    if isequal(nf2,0), s = [-s,s]; end%if
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LOCAL FUNCTION: pscales

function [pr,pm] = pscales(p,nf,ns)
% summary: Parameter perturbation scales

    assert(size(p,2) >= 2,' emgr: min and max parameter required!');

    [pmin,pmax] = bounds(p,2);

    switch nf
        case 1     % Linear centering and scaling
            pr = 0.5 * (pmax + pmin);
            pm = (pmax - pmin) * linspace(0,1.0,ns) + (pmin - pr);

        case 2     % Logarithmic centering and scaling
            lmin = log(pmin);
            lmax = log(pmax);
            pr = real(exp(0.5 * (lmax + lmin)));
            pm = real(exp((lmax - lmin) * linspace(0,1.0,ns) + lmin)) - pr;

        otherwise  % No centering and linear scaling
            pr = pmin;
            pm = (pmax - pmin) * linspace(1.0 / ns,1.0,ns);
    end%switch
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LOCAL FUNCTION: id

function x = id(x,u,p,t)
% summary: (Output) identity functional

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LOCAL FUNCTION: ainv

function x = ainv(m)
% summary: Quadratic complexity approximate inverse matrix

    % Based on truncated Neumann series: X = D^-1 - D^-1 (M - D) D^-1
    D = diag(m);
    k = find(abs(D) > sqrt(eps));
    D(k) = 1.0 ./ D(k);
    x = m .* (-D);
    x = x .* (D');
    x(1:numel(D) + 1:end) = D;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LOCAL FUNCTION: ssp2

function y = ssp2(f,g,t,x0,u,p)
% summary: Low-Storage Strong-Stability-Preserving Second-Order Runge-Kutta

    global STAGES;  % Configurable number of stages for enhanced stability
    if not(isscalar(STAGES)), STAGES = 3; end%if

    dt = t(1);
    nt = floor(t(2) / dt) + 1;

    y0 = g(x0,u(0),p,0);
    y = zeros(numel(y0),nt);  % Pre-allocate trajectory
    y(:,1) = y0;

    xk1 = x0;
    xk2 = x0;
    for k = 2:nt
        tk = (k - 1.5) * dt;
        uk = u(tk);
        for s = 1:(STAGES - 1)
            xk1 = xk1 + (dt / (STAGES - 1)) * f(xk1,uk,p,tk);
        end%for
        xk2 = xk2 + dt * f(xk1,uk,p,tk);
        xk2 = xk2 / STAGES;
        xk2 = xk2 + xk1 * ((STAGES - 1) / STAGES);
        xk1 = xk2;
        y(:,k) = g(xk1,uk,p,tk);
    end%for
end

