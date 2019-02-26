function curios(sys,task,method,options)
%%% summary: curios - Clearing Up Reducibility of Input-Output Systems
%%% project: emgr - EMpirical GRamian Framework ( https://gramian.de )
%%% authors: Christian Himpe ( 0000-0003-2194-6754 )
%%% license: BSD-2-Clause (2019)

    global ODE;
    global STAGES;

    persistent EV;
    persistent ev;
    persistent SV;

%% State-Space Structure

    if(isobject(sys))

        sys = struct(sys);

        if(isfield(sys,'A')), [sys.a] = sys.A; end;
        if(isfield(sys,'B')), [sys.b] = sys.B; end;
        if(isfield(sys,'C')), [sys.c] = sys.C; end;
        if(isfield(sys,'D')), [sys.d] = sys.D; end;
        if(isfield(sys,'E')), [sys.e] = sys.E; end;

        sys.M = size(sys.b,2);
        sys.N = size(sys.a,2);
        sys.Q = size(sys.c,1);
        sys.dt = 1.0/sys.N;
        sys.Tf = 1.0;

        if(isfield(sys,'e') && not(isempty(sys.e)))
            [L,U,P] = lu(sys.e);
            sys.f = @(x,u,p,t) L\(U\(P*(sys.a*x + sys.b*u)));
            sys.F = @(x,u,p,t) L\(U\(P*(sys.a*x))) + sys.c'*u;
        else
            sys.f = @(x,u,p,t) sys.a*x + sys.b*u;
            sys.F = @(x,u,p,t) sys.a'*x + sys.c'*u;
        end

        if(not(isfield(sys,'d')) || isempty(sys.d) || sys.d==0)
            sys.g = @(x,u,p,t) sys.c*x;
        else
            sys.g = @(x,u,p,t) sys.c*x + sys.d*u;
        end
    end

%% Default Arguments

    in = ones(sys.M,1);

    if(not(isfield(sys,'g'))),  sys.g = 1; end
    if(not(isfield(sys,'p'))),  sys.p = 0; end
    if(not(isfield(sys,'q'))),  sys.q = 0; end
    if(not(isfield(sys,'ut'))), sys.ut = 'i'; end
    if(not(isfield(sys,'vt'))), sys.vt = @(t) in*((t<=sys.dt)./sys.dt); end
    if(not(isfield(sys,'us'))), sys.us = 0; end
    if(not(isfield(sys,'xs'))), sys.xs = zeros(sys.N,1); end
    if(not(isfield(sys,'um'))), sys.um = 1; end
    if(not(isfield(sys,'xm'))), sys.xm = 1; end

%% Argument Check

    task = lower(task);
    method = lower(method);

    if(nargin<4 || isempty(options)), options = {'none'}; end

%% Internal Argument setup

    sysdim = [sys.M,sys.N,sys.Q];
    tdisc = [sys.dt,sys.Tf];
    config = zeros(1,12);
    proj = '';
    dp = @mtimes;

%% Utility Library

    picked = @(name) any(strcmp(options,name));
    gtimes = @(m) m*m';
    rcumsum = @(v) flipud(cumsum(flipud(v(:))));

%% Input Library:

    if(picked('step')),   sys.ut = 's'; end
    if(picked('chirp')),  sys.ut = 'c'; end
    if(picked('random')), sys.ut = 'r'; end

%% Kernel Library:

    DEG = 2;

    % Diagonal-ony Pseudo-kernel
    o_diag = @(x,y) sum(x.*y',2);

    % Trace-only pseudo kernel
    o_trac = @(x,y) sum(sum(x.*y'));

    % Time-weighted kernel
    if(picked('tweighted')),  dp = @(x,y) bsxfun(@times,x,(sys.dt * [1:1:size(x,2)]).^DEG) * y; end

    % Polynomial kernel
    if(picked('polynomial')), dp = @(x,y) (x*y).^DEG + 1.0; end

    % Sigmoid kernel
    if(picked('sigmoid')),    dp = @(x,y) tanh(x*y) + 1.0; end

    % Gaussian kernel
    if(picked('gauss')),      dp = @(x,y) exp(-gtimes(x - y')); end

    % Second-order position kernel
    if(picked('position')),   dp = @(x,y) x(1:sys.N/2,:)*y(:,1:sys.N/2); proj = 'secondo'; end

    % Second-order velocity kernel
    if(picked('velocity')),   dp = @(x,y) x(sys.N/2+1:end,:)*y(:,sys.N/2+1:end); proj = 'secondo'; end

%% Solver Library

    if(picked('rk1e')),       STAGES = 1; end
    if(picked('rk45')),       ODE = @rk45e; end

%% Configuration Library

    if(picked('jacobi')),   config(6) = 1; end
    if(picked('scaled')),   config(6) = 2; end
    if(picked('nonsym')),   config(7) = 1; end
    if(picked('active')),   config(8) = 1; end
    if(picked('linpar')),   config(9) = 1; end
    if(picked('logpar')),   config(9) = 2; end
    if(picked('coarse')),   config(10) = 1; end

    if(picked('steady')),   config(1) = 1; end
    if(picked('final')),    config(1) = 2; end
    if(picked('mean')),     config(1) = 3; end
    if(picked('rms')),      config(1) = 4; end
    if(picked('midr')),     config(1) = 5; end
    if(picked('gmean')),     config(1) = 6; end

    if(picked('linear')),   config(2:3) = 1; end
    if(picked('log')),      config(2:3) = 2; end
    if(picked('geom')),     config(2:3) = 3; end
    if(picked('sparse')),   config(2:3) = 4; end

    % Check for subplot configuration
    sc = cellfun(@(c) isnumeric(c) && numel(c)==3,options);
    if(any(sc))
        subpl = options{sc};
        options(sc) = [];
        if(isempty(options)), options = {'none'}; end
    else
        subpl = [];
    end

    if(picked('hold')), holding = 1; else, holding = 0; end

%% Backend Setup

    [EMGR,mdf] = sel();

%% Print Welcome

    fprintf('\n');
    disp('======== curios - Clearing Up Reducibility of Input-Output Systems ========');
    fprintf(' backend: emgr %.1f%s\n',EMGR('version'),mdf)
    fprintf(' %s: %s \n',task,method);
    fprintf(' options: ');
    for k=1:numel(options)-1, fprintf([options{k},', ']); end
    fprintf([options{end},'\n']);
    fprintf(' system dims: %d input(s), %d states, %d output(s) \n',sys.M,sys.N,sys.Q); 

%% Main body

    tic;

    switch(task)

%% State Reduction

        case 'state-reduction'

            switch(method)

                case 'controllability-truncation'

                    WC = EMGR(sys.f,sys.g,sysdim,tdisc,'c',sys.p,config,sys.ut,sys.us,sys.xs,sys.um,sys.xm,dp);
                    [UX,DX,VX] = svd(WC);

                case 'observability-truncation'

                    WO = EMGR(sys.f,sys.g,sysdim,tdisc,'o',sys.p,config,sys.ut,sys.us,sys.xs,sys.um,sys.xm,dp);
                    [UX,DX,VX] = svd(WO);

                case 'linear-direct-truncation'

                    WX = EMGR(sys.f,sys.F,sysdim,tdisc,'y',sys.p,config,sys.ut,sys.us,sys.xs,sys.um,sys.xm,dp);
                    [UX,DX,VX] = balco(WX);

                case 'nonlinear-direct-truncation'

                    if(picked('partitioned'))
                        config(11) = ceil(sys.N/8);
                        K = ceil(sys.N/config(11));
                        for k=1:K
                            config(12) = k;
                            wx{k} = EMGR(sys.f,sys.g,sysdim,tdisc,'x',sys.p,config,sys.ut,sys.us,sys.xs,sys.um,sys.xm,dp);
                        end
                        WX = cell2mat(wx);
                    else
                        WX = EMGR(sys.f,sys.g,sysdim,tdisc,'x',sys.p,config,sys.ut,sys.us,sys.xs,sys.um,sys.xm,dp);
                    end
                    [UX,DX,VX] = balco(WX);

                case 'linear-balanced-truncation'

                    WC = EMGR(sys.f,1,sysdim,tdisc,'c',sys.p,config,sys.ut,sys.us,sys.xs,sys.um,sys.xm,dp);
                    WO = EMGR(sys.F,1,sysdim,tdisc,'c',sys.p,config,sys.ut,sys.us,sys.xs,sys.um,sys.xm,dp);
                    [UX,DX,VX] = balco(WC,WO);

                case 'nonlinear-balanced-truncation'

                    WC = EMGR(sys.f,sys.g,sysdim,tdisc,'c',sys.p,config,sys.ut,sys.us,sys.xs,sys.um,sys.xm,dp);
                    WO = EMGR(sys.f,sys.g,sysdim,tdisc,'o',sys.p,config,sys.ut,sys.us,sys.xs,sys.um,sys.xm,dp);
                    [UX,DX,VX] = balco(WC,WO);

                case 'linear-dominant-subspaces'
                    WX = EMGR(sys.f,sys.F,sysdim,tdisc,'y',sys.p,config,sys.ut,sys.us,sys.xs,sys.um,sys.xm,dp);
                    [UX,DX,VX] = svd(WX);
                    UX = UX*DX;
                    VX = VX*DX;
                    proj = 'dominant';

                case 'nonlinear-dominant-subspaces'
                    WX = EMGR(sys.f,sys.g,sysdim,tdisc,'x',sys.p,config,sys.ut,sys.us,sys.xs,sys.um,sys.xm,dp);
                    [UX,DX,VX] = svd(WX);
                    UX = UX*DX;
                    VX = VX*DX;
                    proj = 'dominant';

                otherwise
                    error('curios: unknown state-reduction method.');
            end

            if(picked('gains')), [UX,DX,VX] = gains(sys,UX,DX,VX); end

            fprintf(' offline time: %.2f s\n',toc);
            [ord,err,nam] = assess(sys,UX,VX,[],[],proj);
            plot_error_1d(ord,err,nam,method,subpl);

            if(not(picked('noscore'))), morscore(ord{1},err{3},sys.N); end;

%% Parameter Reduction

        case 'parameter-reduction'

            switch(method)

                case 'controllability-based'

                    WI = EMGR(sys.f,sys.g,sysdim,tdisc,'s',sys.p,config,sys.ut,sys.us,sys.xs,sys.um,sys.xm,dp);

                case 'observability-based'

                    WI = EMGR(sys.f,sys.g,sysdim,tdisc,'i',sys.p,config,sys.ut,sys.us,sys.xs,sys.um,sys.xm,dp);

                case 'minimality-based'

                    if(picked('partitioned'))
                        config(11) = ceil((sys.N+size(sys.p,1))/8);
                        K = ceil((sys.N+size(sys.p,1))/config(11));
                        for k=1:K
                            config(12) = k;
                            wi{1,k} = EMGR(sys.f,sys.g,sysdim,tdisc,'j',sys.p,config,sys.ut,sys.us,sys.xs,sys.um,sys.xm,dp);
                        end
                        WI{1} = cell2mat(wi(1:ceil(sys.N/config(11))));
                        WI{2} = cell2mat(wi(ceil(sys.N/config(11))+1:end));
                        if(config(10))
                            WI{2} = -0.5*WI{2}'*WI{2};
                        else
                            WI{2} = -0.5*WI{2}'*ainv(WI{1}+WI{1}')*WI{2};
                        end
                    else
                        WI = EMGR(sys.f,sys.g,sysdim,tdisc,'j',sys.p,config,sys.ut,sys.us,sys.xs,sys.um,sys.xm,dp);
                    end

                otherwise
                    error('curios: unknown parameter-reduction method.');
            end

            [UP,DP,VP] = svd(WI{2});
            fprintf(' offline time: %.2f s\n',toc);
            [ord,err,nam] = assess(sys,[],[],UP,UP,proj);
            plot_error_1d(ord,err,nam,method,subpl);

            if(not(picked('noscore'))), morscore(ord{2},err{3},size(sys.p,1)); end;

%% Combined State and Parameter Reduction

        case 'combined-reduction'

            switch(method)

                case 'controllability-based'

                    WI = EMGR(sys.f,sys.g,sysdim,tdisc,'s',sys.p,config,sys.ut,sys.us,sys.xs,sys.um,sys.xm,dp);
                    WO = EMGR(sys.F,sys.g,sysdim,tdisc,'c',sys.p,config,sys.ut,sys.us,sys.xs,sys.um,sys.xm,dp);
                    [UX,DX,VX] = balco(WI{1},WO);

                case 'observability-based'

                    WC = EMGR(sys.f,sys.g,sysdim,tdisc,'c',sys.p,config,sys.ut,sys.us,sys.xs,sys.um,sys.xm,dp);
                    WI = EMGR(sys.f,sys.g,sysdim,tdisc,'i',sys.p,config,sys.ut,sys.us,sys.xs,sys.um,sys.xm,dp);
                    [UX,DX,VX] = balco(WC,WI{1});

                case 'minimality-based'

                    if(picked('partitioned'))
                        config(11) = ceil((sys.N+size(sys.p,1))/8);
                        K = ceil((sys.N+size(sys.p,1))/config(11));
                        for k=1:K
                            config(12) = k;
                            wi{1,k} = EMGR(sys.f,sys.g,sysdim,tdisc,'j',sys.p,config,sys.ut,sys.us,sys.xs,sys.um,sys.xm,dp);
                        end
                        WI{1} = cell2mat(wi(1:ceil(sys.N/config(11))));
                        WI{2} = cell2mat(wi(ceil(sys.N/config(11))+1:end));
                        if(config(10))
                            WI{2} = -0.5*WI{2}'*WI{2};
                        else
                            WI{2} = -0.5*WI{2}'*ainv(WI{1}+WI{1}')*WI{2};
                        end
                    else
                        WI = EMGR(sys.f,sys.g,sysdim,tdisc,'j',sys.p,config,sys.ut,sys.us,sys.xs,sys.um,sys.xm,dp);
                    end
                    if(picked('dominant'))
                        [UX,DX,VX] = svd(WI{1});
                        UX = UX*DX;
                        VX = VX*DX;
                        proj = 'dominant';
                    else
                        [UX,DX,VX] = balco(WI{1});
                    end

                otherwise
                    error('curios: unknown combined-reduction method.');
            end

            [UP,DP,VP] = svd(WI{2});
            fprintf(' offline time: %.2f s\n',toc);
            [ord,err,nam] = assess(sys,UX,VX,UP,UP,proj);
            plot_error_2d(ord,err,nam,method,subpl);

%% Sensitivity Analysis

        case 'sensitivity-analysis'

            switch(method)

                case 'input-state-based'

                    WS = EMGR(sys.f,sys.g,sysdim,tdisc,'s',sys.p,config,sys.ut,sys.us,sys.xs,sys.um,sys.xm);

                case 'input-output-based'

                    config(10) = 1;
                    WS = EMGR(sys.f,sys.g,sysdim,tdisc,'s',sys.p,config,sys.ut,sys.us,sys.xs,sys.um,sys.xm);

                otherwise
                    error('curios: unknown sensitivity-analysis method.');
            end

            fprintf(' offline time: %.2f s\n',toc);
            plot_mag_1d(WS{2}./max(WS{2}),method,'sensitivity',holding,subpl);

%% Parameter Identification

        case 'parameter-identification'

            switch(method)

                case 'state-output-based'

                    WI = EMGR(sys.f,sys.g,sysdim,tdisc,'i',sys.p,config,sys.ut,sys.us,sys.xs,sys.um,sys.xm,dp);

                case 'input-output-based'

                    WI = EMGR(sys.f,sys.g,sysdim,tdisc,'j',sys.p,config,sys.ut,sys.us,sys.xs,sys.um,sys.xm,dp);

                otherwise
                    error('curios: unknown parameter-identification method.');
            end

            EO = sort(abs(eig(WI{2})),'descend');
            fprintf(' offline time: %.2f s\n',toc);
            plot_mag_1d(abs(EO)./max(EO),method,'identifiability',holding,subpl);

%% Decentralized Control

        case 'decentralized-control'

            PM = zeros(sys.M,sys.Q);
            for m = 1:sys.M
                u0 = sparse(m,1,1,sys.M,1);
                f = @(x,u,p,t) sys.f(x,u0*u,p,t);
                for q = 1:sys.Q
                    y0 = sparse(1,q,1,1,sys.Q);

                    switch(method)

                        case 'linear'
                        F = @(x,u,p,t) sys.F(x,y0'*u,p,t);
                        WXij = EMGR(f,F,[1,sys.N,1],tdisc,'y',sys.p,config,1,0,sys.xs,1,sys.xm,o_diag);

                        case 'nonlinear'
                        g = @(x,u,p,t) y0*sys.g(x,u0*u,p,t);
                        WXij = EMGR(f,g,[1,sys.N,1],tdisc,'x',sys.p,config,1,0,sys.xs,1,sys.xm,o_diag);

                        otherwise
                            error('curios: unknown decentralized-control method.');
                    end

                    PM(m,q) = sum(WXij.^2);
                    if(picked('coherence')), PM(m,q) = sum(WXij)^2/PM(m,q); end
                end
            end

            PM = PM - min(PM(:));
            PM = PM./max(PM(:));
            fprintf(' offline time: %.2f s\n',toc);
            plot_iomat(PM,sys.M,sys.Q,method,subpl);

%% Nonlinearity Quantification

        case 'nonlinearity-quantification'

            K = 20;
            nl = zeros(1,K);
            sc = linspace(0.1,2.0,K);

            switch(method)

                case 'input-based'

                    for k = 1:K
                        nl(k) = EMGR(sys.f,sys.g,sysdim,tdisc,'c',sys.p,config,sys.ut,sys.us,sys.xs,sc(k)*sys.um,sys.xm,o_trac);
                    end

                case 'state-based'

                    for k = 1:K
                        nl(k) = EMGR(sys.f,sys.g,sysdim,tdisc,'x',sys.p,config,sys.ut,sys.us,sys.xs,sc(k)*sys.um,sc(k)*sys.xm,o_trac);
                    end

                case 'output-based'

                    for k = 1:K
                        nl(k) = EMGR(sys.f,sys.g,sysdim,tdisc,'o',sys.p,config,sys.ut,sys.us,sys.xs,sys.um,sc(k)*sys.xm,o_trac);
                    end

                otherwise
                    error('curios: unknown nonlinearity-quantification method.');
            end

            nl = abs(nl);
            nl = nl - min(nl);
            nl = nl./max(nl);
            fprintf(' offline time: %.2f s\n',toc);
            plot_mag_1d(nl,method,'nonlinearity',holding,subpl);

%% State Index

        case 'state-index'

            switch(method)

                case 'controllability'

                    w = EMGR(sys.f,sys.g,sysdim,tdisc,'c',sys.p,config,sys.ut,sys.us,sys.xs,sys.um,sys.xm,o_diag);

                case 'observability'

                    w = EMGR(sys.f,sys.g,sysdim,tdisc,'o',sys.p,config,sys.ut,sys.us,sys.xs,sys.um,sys.xm,o_diag);

                case 'minimality'

                    w = EMGR(sys.f,sys.g,sysdim,tdisc,'x',sys.p,config,sys.ut,sys.us,sys.xs,sys.um,sys.xm,o_diag);

                otherwise
                    error('curios: unknown state-index method.');
            end

            fprintf(' offline time: %.2f s\n',toc);
            plot_mag_1d(abs(w)./max(abs(w)),method,task,holding,subpl);

%% System Indices

        case 'system-index'

            if(not(picked('cached')) || isempty(EV))
                WX = EMGR(sys.f,sys.g,sysdim,tdisc,'x',sys.p,config,sys.ut,sys.us,sys.xs,sys.um,sys.xm,dp);
                EV = eig(WX);
                [~,IX] = sort(abs(EV),'descend');
                EV = real(EV(IX));
                ev = EV;
                ev(abs(ev)<(eps*max(abs(ev)))) = 0;
                SV = svd(WX);
            end

            switch(method)

                case 'cauchy-index'

                    in = sum(sign(ev)) - cumsum(sign(ev));

                case 'system-entropy'

                    in = 1.0 - ([1:sys.N]'*log(2*pi*exp(1)) + cumsum(log(abs(EV))))./(sys.N*log(2*pi*exp(1)) + sum(log(abs(EV))));

                case 'system-gain'

                    in = 1.0 - cumsum(EV)./sum(EV);
                    in = in./max(in);

                case 'hinf-bound'

                    in = 2.0 * rcumsum([abs(EV(2:end));0])./sum(abs(EV));
                    in = in./max(in);

                case 'hankel-bound'

                    in = [abs(EV(2:end));0]./sum(abs(EV));
                    in = in./max(in);

                case 'energy-fraction'

                    in = 1.0 - cumsum(SV)./sum(SV);
                    in = in./max(in);

                case 'storage-efficiency'

                    in = cumprod(abs(EV)).^(1.0./[1:numel(EV)]');
                    in = in./max(in);

                case 'ellipsoid-volume'

                    in = sqrt(cumprod(abs(EV)));
                    in = max(1e-16,in./max(in));

                case 'nyquist-area' % = System-Norm

                    in = 1.0 - sqrt(cumsum(EV.^2) ./ sum(EV.^2));
                    in = in./max(in);

                case 'io-coherence'

                    in = 1.0 - abs((cumsum(EV).^2 * sum(EV.^2)) ./ (cumsum(EV.^2) * sum(EV).^2));
                    in = in./max(in);

                otherwise
                    error('curios: unknown system-index method.');
            end

            fprintf(' offline time: %.2f s\n',toc);
            plot_mag_1d(in,method,task,holding,subpl);

        otherwise

            error('Unknown Task!');
    end

%% Clean up

    ODE = [];
    STAGEC = [];

    fprintf('\n');
end

%% LOCAL FUNCTION: select emgr backend
function [EMGR,mdf] = sel()

    if(exist('OCTAVE_VERSION','builtin'))
        assert(exist('emgr_oct')==2,'emgr_oct not found! Get emgr at: https://gramian.de');
        EMGR = @emgr_oct;
        mdf = '-oct';
    elseif(verLessThan('matlab','9.1'))
        assert(exist('emgr_lgc')==2,'emgr_lgc not found! Get emgr at: https://gramian.de');
        EMGR = @emgr_lgc;
        mdf = '-lgc';
    else

        assert(exist('emgr')==2,'emgr not found! Get emgr at: https://gramian.de');
        EMGR = @emgr;
        mdf = '';
    end
end

%% LOCAL FUNCTION: ainv (approximate inverse)
function x = ainv(m)

    d = diag(m);
    k = find(abs(d)>sqrt(eps));
    d(k) = 1.0./d(k);
    x = bsxfun(@times,m,-d);
    x = bsxfun(@times,x,d');
    x(1:numel(d)+1:end) = d;
end

%% LOCAL FUNCTION: balco (Balancing Controllability and Observability)
function [TR,HSV,TL] = balco(WC,WO)

    switch(nargin)

        case 1
            [LO,EX,LC] = eig(WC,'vector');
            [ex,IN] = sort(sqrt(abs(EX)),'descend');
            LO = LO(:,IN)*diag(ex);
            LC = LC(:,IN)*diag(ex);
        case 2
            [LC,EC] = eig(WC,'vector');
            [ec,IN] = sort(sqrt(abs(EC)),'descend');
            [LO,EO] = eig(WO,'vector');
            [eo,IN] = sort(sqrt(abs(EO)),'descend');
            LC = LC(:,IN)*diag(ec);
            LO = LO(:,IN)*diag(eo);
    end

    [U,HSV,V] = svd(LC'*LO);
    HSV = diag(HSV);
    y = sqrt(1.0./(HSV + 1e-14));
    TR = LO*V*diag(y);
    TL = LC*U*diag(y);
end

%% LOCAL FUNCTION: gains (Balanced Gains Sorting)
function [TR,HSV,TL] = gains(sys,tr,hs,tl)

    B = zeros(sys.N,sys.M);
    for k = 1:sys.M
        B(:,k) = sys.f(sys.xs,(1:sys.M==k)',sys.p,0);
    end
    C = sys.g(1,sys.us,sys.p,0);
    HSV = abs(sum((tr'*B).*(C*tl)',2)).*hs;
    [TMP,IX] = sort(HSV,'descend');
    TL = tl(:,IX);
    TR = tr(:,IX);
end

%% LOCAL FUNCTION: assess (Assess Reduced Order Model)
function [orders,errors,names] = assess(sys,UX,VX,UP,VP,m)

    global ODE;

    P = size(sys.p,1);
    Q = size(sys.q,2);

    if(not(isempty(UX)) && not(isempty(UP)))
        if(sys.N>10), xskip = round(sys.N./10); else, xskip = 1; end
        if(P>10), pskip = round(P./10); else, pskip = 1; end
    else
        if(sys.N>100), xskip = round(sys.N./100); else, xskip = 1; end
        if(P>100), pskip = round(P./100); else, pskip = 1; end
    end

    xtest = [1,xskip:xskip:min(size(UX,2),sys.N)-1];
    ptest = [1,pskip:pskip:min(size(UP,2),P)-1];

    if(isempty(UX) || isempty(VX))
        UX = 1;
        VX = 1;
        xtest = 1;
    end

    if(isempty(UP) || isempty(VP))
        UP = 1;
        VP = 1;
        ptest = 1;
    end

    l0norm = @(y) sum(prod(abs(y),1).^(1.0/size(y,1)),2);
    l1norm = @(y) norm(y(:),1);
    l2norm = @(y) norm(y(:),2);
    l8norm = @(y) norm(y(:),Inf);

    l0 = zeros(numel(xtest),numel(ptest));
    l1 = zeros(numel(xtest),numel(ptest));
    l2 = zeros(numel(xtest),numel(ptest));
    l8 = zeros(numel(xtest),numel(ptest));

    for q = 1:Q

        Y = ODE(sys.f,sys.g,[sys.dt,sys.Tf],sys.xs,@(t) sys.us + sys.vt(t),sys.q(:,q));
        n0 = l0norm(Y);
        n1 = l1norm(Y);
        n2 = l2norm(Y);
        n8 = l8norm(Y);

        for n = xtest

            switch(m)

                case 'dominant'
                    [ux,~,~] = svd([UX(:,1:n),VX(:,1:n)],'econ');
                    ux = ux(:,1:n);
                    vx = ux';

                case 'secondo'
                    nn = numel(sys.xs)/2;
                    ux = [UX(:,1:n),zeros(nn,n);zeros(nn,n),UX(:,1:n)];
                    vx = [VX(:,1:n),zeros(nn,n);zeros(nn,n),VX(:,1:n)]';

                otherwise
                    ux = UX(:,1:n);
                    vx = VX(:,1:n)';
            end

            ix = find(xtest==n);

            for p = ptest

                up = UP(:,1:p);
                vp = VP(:,1:p)';

                ip = find(ptest==p);

                y = ODE(@(x,u,p,t) vx*sys.f(ux*x,u,up*p,t), ...
                        @(x,u,p,t)    sys.g(ux*x,u,up*p,t), ...
                        [sys.dt,sys.Tf], vx*sys.xs, @(t) sys.us + sys.vt(t), vp*sys.q(:,q));

                e0 = l0norm(Y - y) / n0;
                e1 = l1norm(Y - y) / n1;
                e2 = l2norm(Y - y) / n2;
                e8 = l8norm(Y - y) / n8;

                if(Q > 1)
                    e0 = e0 .* e0;
                    e1 = e1 .* e1;
                    e2 = e2 .* e2;
                    e8 = e8 .* e8;
                end

                l0(ix,ip) = l0(ix,ip) + e0;
                l1(ix,ip) = l1(ix,ip) + e1;
                l2(ix,ip) = l2(ix,ip) + e2;
                l8(ix,ip) = l8(ix,ip) + e8;
            end
        end
    end

    if(Q>1)
        l0 = sqrt(l0 ./ Q);
        l1 = sqrt(l1 ./ Q);
        l2 = sqrt(l2 ./ Q);
        l8 = sqrt(l8 ./ Q);
    end

    l0 = min(l0,1.0);
    l1 = min(l1,1.0);
    l2 = min(l2,1.0);
    l8 = min(l8,1.0);

    orders = {xtest, ptest};
    errors = {l0, l1, l2, l8};
    names = {'L_1 Error', 'L_2 Error', 'L_\infty Error', 'L_0 Error'};
end

%% LOCAL FUNCTION: morscore (Model Reduction Scoring)
function morscore(orders,errors,N)

    nx = orders./(N-1);
    ny = log10(errors)./16.0 + 1.0;
    if(nx(end)~=1), nx(end+1) = 1; ny(end+1) = ny(end); end;
    ms = 1.0 - trapz(nx(:),ny(:));
    fprintf(' MOR score: %.2f',ms);
    text(0.5,0.5,num2str(ms,'%.2f'),'Units','normalized');
end

%% LOCAL FUNCTION: rk45e (Embedded Runge Kutta 4th / 5th Solver)
function y = rk45e(f,g,t,x0,u,p)

    [S,x] = ode45(@(t,x) f(x,u(t),p,t),[0,t(2)],x0,'InitialStep',t(1)); % Compute State Trajectory

    K = numel(S);
    z = g(x(1,:)',u(S(1)),p,S(1));
    z(end,K) = 0;

    for k = 2:K
        tk = S(k);
        z(:,k) = g(x(k,:)',u(tk),p,tk); % Compute Output Trajectory
    end

    y = interp1(S,z',0:t(1):t(2))';
end

%% LOCAL FUNCTION: plot_error_1d (1D Error Plot)
function plot_error_1d(orders,errors,names,ident,subpl)

    if(orders{1}==1), dom = orders{2}; else, dom = orders{1}; end

    if(isempty(subpl))
        figure('Name',ident,'NumberTitle','off');
    else
        if(subpl(3)==1), figure; end
        subplot(subpl(1),subpl(2),subpl(3));
    end

    semilogy(dom,errors{2},'r','linewidth',2); hold on;
    semilogy(dom,errors{3},'g','linewidth',2);
    semilogy(dom,errors{4},'b','linewidth',2);
    semilogy(dom,errors{1},'k--','linewidth',2); hold off;
    xlim([dom(1),dom(end)]);
    ylim([1e-16,1]);
    pbaspect([2,1,1]);

    if(isempty(subpl))
        xlabel('Reduced Dimension');
        ylabel('Relative Error');
        legend(names{1},names{2},names{3},names{4},'location','northeast');
    else
        set(gca,'Ytick',[]);
        set(gca,'Xtick',[]);
        if(subpl(3)<=subpl(2)), title(ident); end
        set([gca; findall(gca,'Type','text')],'FontSize',4);
    end
    set(gca,'YGrid','on');
end

%% LOCAL FUNCTION: plot_error2d (2D Error Plot)
function plot_error_2d(orders,errors,names,ident,subpl)

    if(isempty(subpl))
        figure('Name',ident,'NumberTitle','off');
    else
        if(subpl(3)==1), figure; end
        subplot(subpl(1),subpl(2),subpl(3));
    end

    h = surf(orders{1},orders{2},min(1.0,errors{3}));
    set(gca,'ZScale','log');
    set(h,'CData',log10(get(h,'CData')));
    set(gca,'CLim',log10(get(gca,'ZLim')));
    view(135,15);

    if(isempty(subpl))
        zl = zlim();
        zlim([zl(1),1]);
        ylabel('State Dimension')
        xlabel('Parameter Dimension');
        zlabel('Relative Error');
    else
        zlim([1e-16,1]);
        set(gca,'Ytick',[]);
        set(gca,'Xtick',[]);
        if(subpl(3)<=subpl(2)), title(ident); end
        set([gca; findall(gca,'Type','text')],'FontSize',4);
    end
end

%% LOCAL FUNCTION: plot_iomat (Plot Participation Matrix)
function plot_iomat(PM,M,Q,ident,subpl)

    if(isempty(subpl))
        figure('Name',ident,'NumberTitle','off');
    else
        if(subpl(3)==1), figure; end
        subplot(subpl(1),subpl(2),subpl(3));
    end

    imagesc(PM);

    if(isempty(subpl))
        colorbar;
        set(gca,'XTick',1:1:M,'YTick',1:1:Q);
    else
        set(gca,'Ytick',[]);
        set(gca,'Xtick',[]);
        set(gca,'Ztick',[]);
        title(ident);
    end
end

%% LOCAL FUNCTION: plot_mag_1d (1D Magnitude Plot)
function plot_mag_1d(mag,name,ident,holding,subpl)

    persistent c;

    if(isempty(c) || (not(holding) && not(isempty(c))) )
        c = 1;
    end

    dom = numel(mag);

    if(holding)
        col = lines(11);
        hold on;
    else
        col = zeros(11,3);

        if(isempty(subpl))
            figure('Name',ident,'NumberTitle','off');
        else
            if(subpl(3)==1), figure; end
            subplot(subpl(1),subpl(2),subpl(3));
        end
    end

    if(strcmp(name,'cauchy-index') || strcmp(name,'system-entropy'))

        plot(1:dom,mag,'Color',col(c,:),'linewidth',2,'DisplayName',name);
    else
        mag(mag==0) = min(mag(mag>0));
        semilogy(1:dom,mag,'Color',col(c,:),'linewidth',2,'DisplayName',name);
    end

    if(holding)
        hold off;
        c = c + 1;
    end

    if(isempty(subpl))
        legend show
        set(legend,'location','northeast');
    else
        title(ident);
    end

    xlim([1,dom]);
    yl = ylim();
    ylim([min([1e-4;yl(1);mag(:)]),max([1;yl(2);mag(:)])]);
    pbaspect([2,1,1]);
    set(gca,'YGrid','on');
end
