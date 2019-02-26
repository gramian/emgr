function W = emgr(f,g,s,t,w,pr,nf,ut,us,xs,um,xm,dp)
%% emgr - EMpirical GRamian Framework
%
%  project: emgr ( https://gramian.de )
%  version: 5.7 ( 2019-02-26 )
%  authors: Christian Himpe ( 0000-0003-2194-6754 )
%  license: BSD-2-Clause License ( opensource.org/licenses/BSD-2-Clause )
%  summary: Empirical system Gramians for (nonlinear) input-output systems.
%
% USAGE:
%
%  W = emgr(f,g,s,t,w,[pr],[nf],[ut],[us],[xs],[um],[xm],[dp])
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
% ARGUMENTS:
%
%   f {handle} vector field handle: x' = f(x,u,p,t)
%   g {handle} output function handle: y = g(x,u,p,t)
%   s {vector} system dimensions: [inputs,states,outputs]
%   t {vector} time discretization: [time-step,time-horizon]
%   w  {char}  single character encoding gramian type:
%    * 'c' empirical controllability gramian (Wc)
%    * 'o' empirical observability gramian (Wo)
%    * 'x' empirical cross gramian (Wx aka Wco or Xcg)
%    * 'y' empirical linear cross gramian (Wy)
%    * 's' empirical sensitivity gramian (Ws)
%    * 'i' empirical identifiability gramian (Wi)
%    * 'j' empirical joint gramian (Wj)
%  pr {matrix|0} parameters, each column is one set
%  nf {vector|0} option flags, twelve component vector, default zero:
%    * center: none(0), steady(1), last(2), mean(3), rms(4), midr(5), geom(6)
%    * input scales: single(0), linear(1), geom(2), log(3), sparse(4)
%    * state scales: single(0), linear(1), geom(2), log(3), sparse(4)
%    * input rotations: unit(0), single(1)
%    * state rotations: unit(0), single(1)
%    * normalization (only: Wc, Wo, Wx, Wy): none(0), Jacobi(1), steady(2)
%    * state gramian variant:
%      * cross gramian type (only: Wx, Wy, Wj): regular(0), non-symmetric(1)
%      * observability gramian type (only: Wo, Wi): regular(0), averaged(1)
%    * extra input (only: Wo, Wx, Ws, Wi, Wj): none(0), yes(1)
%    * parameter centering (only: Ws, Wi, Wj): none(0), linear(1), log(2)
%    * parameter gramian variant:
%      * averaging type (only: Ws): input-state(0), input-output(1)
%      * Schur-complement (only: Wi, Wj): detailed(0), approximate(1)
%    * cross gramian partition size (only: Wx, Wj): full(0), partitioned(<N)
%    * cross gramian partition index (only: Wx, Wj): partition(>0)
%  ut {handle|'i'} input function handle: u_t = ut(t) or character:
%    * 'i' delta impulse input
%    * 's' step input / load vector / source term
%    * 'c' decaying exponential chirp input
%    * 'r' pseudo-random binary input
%  us {vector|0} steady-state input
%  xs {vector|0} steady-state and nominal initial state x_0
%  um {matrix|1} input scales
%  xm {matrix|1} initial-state scales
%  dp {handle|@mtimes} inner product handle: xy = dp(x,y)
%
% RETURNS:
%
%  W {matrix} Gramian Matrix (for: Wc, Wo, Wx, Wy)
%  W  {cell}  [State-, Parameter-] Gramian (for: Ws, Wi, Wj)
%
% CITE AS:
%
%  C. Himpe (2019). emgr - EMpirical GRamian Framework (Version 5.7)
%  [Software]. Available from https://gramian.de . doi:10.5281/zenodo.2577980
%
% KEYWORDS:
%
%  model reduction, system gramians, empirical gramians, cross gramian, MOR
%
% SEE ALSO: gram (Control System Toolbox)

    % Integrator Handle
    global ODE;
    if not(isa(ODE,'function_handle')), ODE = @ssp2; end%if

    % Version Info
    if strcmp(f,'version'), W = 5.7; return; end%if

    % Default Arguments
    if (nargin <  6) || isempty(pr), pr = 0.0; end%if
    if (nargin <  7) || isempty(nf), nf = 0.0; end%if
    if (nargin <  8) || isempty(ut), ut = 'i'; end%if
    if (nargin <  9) || isempty(us), us = 0.0; end%if
    if (nargin < 10) || isempty(xs), xs = 0.0; end%if
    if (nargin < 11) || isempty(um), um = 1.0; end%if
    if (nargin < 12) || isempty(xm), xm = 1.0; end%if
    if (nargin < 13) || isempty(dp), dp = @mtimes; end%if

%% SETUP

    % System Dimensions
    M = s(1);                 % Number of inputs
    N = s(2);                 % Number of states
    Q = s(3);                 % Number of outputs
    A = 0;                    % Number of augmented parameter states
    P = size(pr,1);           % Dimension of parameter and number of sets
    K = size(pr,2);           % Number of parameter-sets

    % Time Discretization
    dt = t(1);                % Time-step width
    Tf = t(2);                % Time horizon
    nt = floor(Tf / dt) + 1;  % Number of time-steps plus initial value

    % Lazy Output Functional
    if isnumeric(g) && (g == 1), g = @id; Q = N; end%if

    % Augmented Parameter-States (set by parameter gramians)
    if (numel(s) == 4), A = s(4); end%if

    % Pad Flag Vector
    if numel(nf) < 12, nf(12) = 0; end%if

    % Built-in Input Functions
    if not(isa(ut,'function_handle'))
        switch lower(ut)

            case 's'  % Step Input
                ut = @(t) 1;

            case 'c'  % Decaying Exponential Chirp Input
                a0 = (2.0 * pi) / (4.0 * dt) * Tf / log(4.0 * (dt / Tf));
                b0 = (4.0 * (dt / Tf)) ^ (1.0 / Tf);
                ut = @(t) 0.5 * cos(a0 * (b0 ^ t - 1.0)) + 0.5;

            case 'r'  % Pseudo-Random Binary Input
                ut = @(t) randi([0,1],1,1);

            otherwise % Delta Impulse Input
                ut = @(t) (t <= dt) / dt;
        end%switch
    end%if

    % Lazy Optional Arguments
    if isscalar(us), us = repmat(us,M,1); end%if
    if isscalar(xs), xs = repmat(xs,N,1); end%if
    if isscalar(um), um = repmat(um,M,1); end%if
    if isscalar(xm), xm = repmat(xm,N,1); end%if

    % Gramian Normalization
    if nf(6) && (A == 0)
        switch nf(6)

            case 1 % Jacobi-type preconditioner
                NF = nf; NF(6) = 0;
                DP = @(x,y) sum(x .* y',2); % Diagonal-only pseudo-kernel
                WT = emgr(f,g,s,t,w,pr,NF,ut,us,xs,um,xm,DP);
                TX = sqrt(abs(WT));

            case 2 % Steady-state preconditioner
                TX = xs;
        end%switch
        TX(abs(TX) < sqrt(eps)) = 1;
        f = @(x,u,p,t) f(TX .* x,u,p,t) ./ TX;
        g = @(x,u,p,t) g(TX .* x,u,p,t);
        xs = xs ./ TX;
    end%if

    % Non-symmetric cross Gramian or average observability Gramian
    if nf(7), R = 1; else, R = Q; end%if

    % Extra Input
    if nf(8), up = @(t) us + ut(t); else, up = @(t) us; end%if

    % Scale Sampling
    if size(um,2) == 1, um = um * scales(nf(2),nf(4)); end%if
    if size(xm,2) == 1, vm = xm(1:Q) * scales(nf(2),nf(4)); end%if
    if size(xm,2) == 1, xm = xm * scales(nf(3),nf(5)); end%if

    C = size(um,2); % Number of input scales sets
    D = size(xm,2); % Number of state scales sets

%% GRAMIAN COMPUTATION

    switch lower(w) % Empirical system gramian types

        % Common Layout:
        %   For each {parameter set, scale, input/state/parameter component}:
        %     Perturb, simulate, center, normalize, accumulate
        %   Assemble, normalize, post-process
        %   Parameter gramians call state gramians

        case 'c' % Empirical Controllability Gramian

            W = 0; % Reserve gramian variable
            for k = 1:K
                for c = 1:C
                    for m = find(um(:,c))' % parfor
                        em = sparse(m + M * (A > 0),1,um(m,c),M + P,1);
                        umc = @(t) up(t) + ut(t) .* em(1:M);
                        pmc = pr(:,k) + em(M + 1:end);
                        x = ODE(f,@id,t,xs,umc,pmc);
                        x = x - avg(x,nf(1),xs);
                        x = x / um(m,c);
                        if A > 0
                            W = W + em(M + 1:end) * dp(x,x');
                        else
                            W = W + dp(x,x');
                        end%if
                    end%for
                end%for
            end%for
            W = W * (dt / (C * K));

        case 'o' % Empirical Observability Gramian

            W = 0; % Reserve gramian variable
            o = zeros(R * nt,N + A); % Pre-allocate observability matrix
            for k = 1:K
                for d = 1:D
                    for n = find(xm(:,d))' % parfor
                        en = sparse(n,1,xm(n,d),N + P,1);
                        xnd = xs + en(1:N);
                        pnd = pr(:,k) + en(N + 1:end);
                        ys = g(xs,us,pnd,0);
                        y = ODE(f,g,t,xnd,up,pnd);
                        y = y - avg(y,nf(1),ys);
                        y = y / xm(n,d);
                        if nf(7) % Average observability gramian
                            o(:,n) = sum(y,1)';
                        else     % Regular observability gramian
                            o(:,n) = y(:);
                        end%if
                    end%for
                    W = W + dp(o',o);
                end%for
            end%for
            W = W * (dt / (D * K));

        case 'x' % Empirical Cross Gramian

            assert((M == Q) || nf(7),'emgr: non-square system!');

            i0 = 1;
            i1 = N + A;

            % Partitioned cross gramian
            if nf(11) > 0
                sp = round(nf(11));      % Partition size
                ip = round(nf(12));      % Partition index
                i0 = i0 + (ip - 1) * sp; % Start index
                i1 = min(i0 + sp - 1,N); % End index
                if i0 > N
                    i0 = i0 - (ceil(N / sp) * sp - N);
                    i1 = min(i0 + sp - 1,N + A);
                end%if

                if (ip < 1) || (i0 > i1) || (i0 < 0), W = 0; return; end%if
            end%if

            W = 0; % Reserve gramian variable
            o = zeros(nt,i1 - i0 + 1,R); % Pre-allocate observability 3-tensor
            for k = 1:K
                for d = 1:D
                    for n = find(xm(i0:i1,d))'
                        en = sparse(i0 - 1 + n,1,xm(i0 - 1 + n,d),N + P,1);
                        xnd = xs + en(1:N);
                        pnd = pr(:,k) + en(N + 1:end);
                        ys = g(xs,us,pnd,0);
                        y = ODE(f,g,t,xnd,up,pnd);
                        y = y - avg(y,nf(1),ys);
                        y = y / xm(i0 - 1 + n,d);
                        if nf(7) % Non-symmetric cross gramian
                            o(:,n,1) = sum(y,1)';
                        else     % Regular cross gramian
                            o(:,n,:) = y';
                        end%if
                    end%for
                    for c = 1:C % parfor
                        for m = find(um(:,c))'
                            em = sparse(m,1,um(m,c),M,1);
                            umc = @(t) us + ut(t) .* em;
                            x = ODE(f,@id,t,xs,umc,pr(:,k));
                            x = x - avg(x,nf(1),xs);
                            x = x / um(m,c);
                            if nf(7) % Non-symmetric cross gramian
                                W = W + dp(x,o(:,:,1));
                            else     % Regular cross gramian
                                W = W + dp(x,o(:,:,m));
                            end%if
                        end%for
                    end%for
                end%for
            end%for
            W = W * (dt / (C * D * K));

        case 'y' % Empirical Linear Cross Gramian

            assert((M == Q) || nf(7),'emgr: non-square system!');
            assert(C == size(vm,2),'emgr: scale count mismatch!');

            W = 0; % Reserve gramian variable
            a = zeros(nt,N,R); % Pre-allocate adjoint 3-tensor
            for k = 1:K
                for c = 1:C
                    for q = find(vm(:,c))'
                        em = sparse(q,1,vm(q,c),Q,1);
                        vqc = @(t) us + ut(t) .* em;
                        z = ODE(g,@id,t,xs,vqc,pr(:,k));
                        z = z - avg(z,nf(1),xs);
                        z = z / vm(q,c);
                        if nf(7) % Non-symmetric cross gramian
                            a(:,:,1) = a(:,:,1) + z';
                        else     % Regular cross gramian
                            a(:,:,q) = z';
                        end%if
                    end%for
                    for m = find(um(:,c))' % parfor
                        em = sparse(m,1,um(m,c),M,1);
                        umc = @(t) us + ut(t) .* em;
                        x = ODE(f,@id,t,xs,umc,pr(:,k));
                        x = x - avg(x,nf(1),xs);
                        x = x / um(m,c);
                        if nf(7) % Non-symmetric cross gramian
                            W = W + dp(x,a(:,:,1));
                        else     % Regular cross gramian
                            W = W + dp(x,a(:,:,m));
                        end%if
                    end%for
                end%for
            end%for
            W = W * (dt / (C * K));

        case 's' % Empirical Sensitivity Gramian

            [pr,pm] = pscales(pr,nf(9),C);
            W{1} = emgr(f,g,[M,N,Q],t,'c',pr,nf,ut,us,xs,um,xm,dp);
            if not(nf(10)) % Input-state sensitivty gramian
                DP = @(x,y) sum(sum(x .* y'));      % Trace pseudo-kernel
            else           % Input-output sensitivity gramian
                DP = @(x,y) sum(reshape(y,Q,[]));   % Custom pseudo-kernel
                Y = emgr(f,g,[M,N,Q],t,'o',pr,nf,ut,us,xs,um,xm,DP);
                DP = @(x,y) abs(sum(y(:) .* Y(:))); % Custom pseudo-kernel
            end%if
            W{2} = emgr(f,g,[M,N,Q,P],t,'c',pr,nf,ut,us,xs,pm,xm,DP);

        case 'i' % Empirical Augmented Observability Gramian

            [pr,pm] = pscales(pr,nf(9),D);
            V = emgr(f,g,[M,N,Q,P],t,'o',pr,nf,ut,us,xs,um,[xm;pm],dp);
            W{1} = V(1:N,1:N);         % Observability gramian
            WM = V(1:N,N + 1:N + P);
            W{2} = V(N + 1:N + P,N + 1:N + P); % Identifiability gramian
            if not(nf(10))
                W{2} = W{2} - (WM' * ainv(W{1}) * WM);
            end%if

        case 'j' % Empirical Joint Gramian

            [pr,pm] = pscales(pr,nf(9),D);
            V = emgr(f,g,[M,N,Q,P],t,'x',pr,nf,ut,us,xs,um,[xm;pm],dp);
            if nf(11), W = V; return; end%if % Joint gramian partition
            W{1} = V(1:N,1:N);               % Cross gramian
            WM = V(1:N,N + 1:N + P);
            if not(nf(10))                   % Cross-identifiability gramian
                W{2} = -0.5 * (WM' * ainv(W{1} + W{1}') * WM);
            else
                W{2} = -0.5 * (WM' * WM);
            end%if

        otherwise

            error('emgr: unknown gramian type!');
    end%switch
end

%% LOCAL FUNCTION: scales
function s = scales(nf1,nf2)
%  summary: Input and initial state perturbation scales

    switch nf1

        case 1 % Linear
            s = [0.25, 0.50, 0.75, 1.0];

        case 2 % Geometric
            s = [0.125, 0.25, 0.5, 1.0];

        case 3 % Logarithmic
            s = [0.001, 0.01, 0.1, 1.0];

        case 4 % Sparse
            s = [0.01, 0.50, 0.99, 1.0];

        otherwise % One
            s = 1.0;
    end%switch

    if nf2 == 0, s = [-s,s]; end%if
end

%% LOCAL FUNCTION: pscales
function [pr,pm] = pscales(p,nf,ns)
%  summary: Parameter perturbation scales

    assert(size(p,2) >= 2,'emgr: min and max parameter required!');

    pmin = min(p,[],2);
    pmax = max(p,[],2);

    switch nf

        case 1 % Linear centering and scales
            pr = 0.5 * (pmax + pmin);
            pm = (pmax - pmin) * linspace(0,1.0,ns) + (pmin - pr);

        case 2 % Logarithmic centering and scales
            lmin = log(pmin);
            lmax = log(pmax);
            pr = real(exp(0.5 * (lmax + lmin)));
            pm = real(exp((lmax - lmin) * linspace(0,1.0,ns) + lmin)) - pr;

        otherwise % No centering and linear scales
            pr = pmin;
            pm = (pmax - pmin) * linspace(1.0 / ns,1.0,ns);
    end%switch
end

%% LOCAL FUNCTION: id
function x = id(x,u,p,t)
%  summary: Output identity function
    ;
end

%% LOCAL FUNCTION: avg
function a = avg(z,nf,zs)
%  summary: State and output trajectory centering

    switch nf

        case 1 % Steady state / output
            a = zs;

        case 2 % Final state / output
            a = z(:,end);

        case 3 % Temporal mean state / output
            a = mean(z,2);

        case 4 % Temporal root-mean-square state / output
            a = sqrt(mean(z .* z,2));

        case 5 % Midrange state / output
            a = 0.5 * (max(z,[],2) - min(z,[],2));

        case 6 % Geometric mean state / output
            a = prod(sign(z),2) .* prod(abs(z),2) .^ (1.0 / size(z,2));

        otherwise % None
            a = 0;
    end%switch
end

%% LOCAL FUNCTION: ainv
function x = ainv(m)
%  summary: Quadratic complexity approximate inverse matrix

    d = diag(m);
    k = find(abs(d) > sqrt(eps));
    d(k) = 1.0 ./ d(k);
    x = m .* (-d);
    x = x .* (d');
    x(1:numel(d) + 1:end) = d;
end

%% LOCAL FUNCTION: ssp2
function y = ssp2(f,g,t,x0,u,p)
%  summary: Low-Storage Strong-Stability-Preserving Second-Order Runge-Kutta

    % Configurable number of stages for enhanced stability
    global STAGES;
    if not(isscalar(STAGES)), STAGES = 3; end%if

    dt = t(1);
    nt = floor(t(2) / dt) + 1;

    y0 = g(x0,u(0),p,0);
    Q = numel(y0);   % Q = N when g = id
    y = zeros(Q,nt); % Pre-allocate trajectory
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
