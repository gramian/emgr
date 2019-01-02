function curios(sys,task,method,options)
%%% summary: curios - Clearing Up Reducibility of Input-Output Systems
%%% project: emgr - EMpirical GRamian Framework ( https://gramian.de )
%%% authors: Christian Himpe ( 0000-0003-2194-6754 )
%%% license: BSD-2-Clause (2019)

    global ODE;

    persistent WX;

%% Default Arguments

    if(not(isfield(sys,'g'))),  sys.g = 1; end
    if(not(isfield(sys,'p'))),  sys.p = 0; end
    if(not(isfield(sys,'q'))),  sys.q = 0; end
    if(not(isfield(sys,'ut'))), sys.ut = 1; end
    if(not(isfield(sys,'vt'))), sys.vt = @(t) ones(sys.M,1)*(t<=sys.h)./sys.h; end
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
    tdisc = [sys.h,sys.T];
    config = zeros(1,12);
    proj = '';
    dp = @mtimes;

%% Utility Library

    picked = @(name) any(strcmp(options,name));
    gtimes = @(m) m*m';
    rcumsum = @(v) flipud(cumsum(flipud(v(:))));

%% Kernel Library:

    DEG = 2;

    % Diagonal-ony Pseudo-kernel
    o_diag = @(x,y) sum(x.*y',2);

    % Trace-only pseudo kernel
    o_trac = @(x,y) sum(sum(x.*y'));

    % Time-weighted kernel
    if(picked('tweighted')),  dp = @(x,y) bsxfun(@times,x,(sys.h * [1:1:size(x,2)]).^DEG) * y; end

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
    if(picked('midrange')), config(1) = 5; end

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

                case 'linear-dominant'
                    WX = EMGR(sys.f,sys.F,sysdim,tdisc,'y',sys.p,config,sys.ut,sys.us,sys.xs,sys.um,sys.xm,dp);
                    [UX,DX,VX] = svd(WX);
                    UX = bsxfun(@times,UX,diag(DX)');
                    VX = bsxfun(@times,VX,diag(DX)');
                    proj = 'dominant';

                case 'nonlinear-dominant'
                    WX = EMGR(sys.f,sys.g,sysdim,tdisc,'x',sys.p,config,sys.ut,sys.us,sys.xs,sys.um,sys.xm,dp);
                    [UX,DX,VX] = svd(WX);
                    UX = bsxfun(@times,UX,diag(DX)');
                    VX = bsxfun(@times,VX,diag(DX)');
                    proj = 'dominant';
            end

            if(picked('gains')), [UX,DX,VX] = gains(sys,UX,DX,VX); end
            fprintf(' offline time: %.2f s\n',toc);
            [ord,err,nam] = assess(sys,UX,VX,[],[],proj);
            plot_error_1d(ord,err,nam,method,subpl);
            %morscore(ord{1},err{3},sys.N); % EXPERIMENTAL

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
                        aelse
                            WI{2} = -0.5*WI{2}'*ainv(WI{1}+WI{1}')*WI{2};
                        end
                    else
                        WI = EMGR(sys.f,sys.g,sysdim,tdisc,'j',sys.p,config,sys.ut,sys.us,sys.xs,sys.um,sys.xm,dp);
                    end
            end

            [UP,DP,VP] = svd(WI{2});
            fprintf(' offline time: %.2f s\n',toc);
            [ord,err,nam] = assess(sys,[],[],UP,UP,proj);
            plot_error_1d(ord,err,nam,method,subpl);
            %morscore(ord{1},err{3},size(sys.p)); % EXPERIMENTAL

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
                        UX = bsxfun(@times,UX,diag(DX)');
                        VX = bsxfun(@times,VX,diag(DX)');
                        proj = 'dominant';
                    else
                        [UX,DX,VX] = balco(WI{1});
                    end
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
            end

            fprintf(' offline time: %.2f s\n',toc);
            plot_mag_1d(WS{2}./sum(WS{2}),method,'sensitivity',holding,subpl);

%% Parameter Identification

        case 'parameter-identification'

            switch(method)

                case 'state-output-based'

                    WI = EMGR(sys.f,sys.g,sysdim,tdisc,'i',sys.p,config,sys.ut,sys.us,sys.xs,sys.um,sys.xm,dp);

                case 'input-output-based'

                    WI = EMGR(sys.f,sys.g,sysdim,tdisc,'j',sys.p,config,sys.ut,sys.us,sys.xs,sys.um,sys.xm,dp);
            end

            EO = sort(abs(eig(WI{2})),'descend');
            fprintf(' offline time: %.2f s\n',toc);
            plot_mag_1d(abs(EO)./sum(EO),method,'identifiability',holding,subpl);

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
            end

            fprintf(' offline time: %.2f s\n',toc);
            plot_mag_1d(abs(w)./sum(abs(w)),method,task,holding,subpl);

%% System Indices

        case 'system-index'

            if(not(picked('cached')) || isempty(WX))
                WX = EMGR(sys.f,sys.g,sysdim,tdisc,'x',sys.p,config,sys.ut,sys.us,sys.xs,sys.um,sys.xm,dp);
            end
            ev = eig(WX);
            [EV,IX] = sort(abs(ev),'descend');
            ev = ev(IX);
            SV = svd(WX);

            switch(method)

                case 'cauchy-index'

                    in = sum(sign(real(ev))) - cumsum(sign(real(ev)));

                case 'system-entropy'

                    in = 1.0 - ([1:sys.N]'*log(2*pi*exp(1)) + cumsum(log(EV)))./(sys.N*log(2*pi*exp(1)) + sum(log(EV)));

                case 'system-gain'

                    in = 1.0 - cumsum(EV)./sum(EV);

                case 'hinf-bound'

                    in = 2.0 * rcumsum([EV(2:end);0])./sum(EV);

                case 'hankel-bound'

                    in = [EV(2:end);0]./sum(EV);

                case 'system-symmetry'

                    in = cumsum(SV.^2)./sum(SV.^2) - cumsum(EV.^2)./sum(EV.^2);

                case 'energy-fraction'

                    in = 1.0 - cumsum(SV)./sum(SV);

                case 'storage-efficiency'

                    in = abs(cumprod(SV).^(1.0./numel(SV)));

                case 'nyquist-area'

                    in = 1.0 - sqrt(pi*cumsum(EV.^2))./sqrt(pi*sum(EV.^2));

                case 'robustness-index'

                    in = abs(2.0 - (4.0/pi) * atan(sqrt(cumsum(EV)./rcumsum(EV))));

                case 'recoverability-index'

                    in = fliplr(EV)./max(EV);

                case 'io-coherence'

                    in = abs(cumsum(ev).^2./cumsum(ev.^2));
                    in = 1.0 - in./max(in);
            end

            fprintf(' offline time: %.2f s\n',toc);
            plot_mag_1d(in,method,task,holding,subpl);

        otherwise

            error('Unknown Task!');
    end

%% Clean up

    ODE = [];

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
        B(:,k) = sys.f(sys.xs,(1:sys.M==1)',sys.p,0);
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

        Y = ODE(sys.f,sys.g,[sys.h,sys.T],sys.xs,@(t) sys.us + sys.vt(t),sys.q(:,q));
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
                        [sys.h,sys.T], vx*sys.xs, @(t) sys.us + sys.vt(t), vp*sys.q(:,q));

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

    orders = {xtest, ptest};
    errors = {l0, l1, l2, l8};
    names = {'L_1 Error', 'L_2 Error', 'L_\infty Error', 'L_0 Error'};
end

%% LOCAL FUNCTION: morscore (MORscore [Experimental; determine weighting])
function morscore(orders,errors,N)

    sigm = @(x,L,k,z) L ./ (1.0 + exp(-k * (x - z)));
    obj = @(p) norm(log10(errors(:))./16.0 - sigm(orders(:),p(1),p(2),p(3)),2);
    depvelshi = fminunc(obj,[-0.5,1,0],optimset('Display','off'));
    score = norm(depvelshi,2);
    fprintf(' MOR score: %.2f',score);
    text(0.5,0.5,num2str(score,'%.2f'),'Units','normalized');
end

%% LOCAL FUNCTION: rk45e (Custom Embedded Runge Kutta 4th / 5th Solver)
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
        set([gca; findall(gca,'Type','text')],'FontSize',6);
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
    zl = zlim();
    zlim([zl(1),1]);
    set(h,'CData',log10(get(h,'CData')));
    set(gca,'CLim',log10(get(gca,'ZLim')));
    view(135,25);

    if(isempty(subpl))
        ylabel('State Dimension')
        xlabel('Parameter Dimension');
        zlabel('Relative Error');
    else
        set(gca,'Ytick',[]);
        set(gca,'Xtick',[]);
        if(subpl(3)<=subpl(2)), title(ident); end
        set([gca; findall(gca,'Type','text')],'FontSize',6);
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
