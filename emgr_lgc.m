function W = emgr_lgc(f,g,s,t,w,pr,nf,ut,us,xs,um,xm,dp)
%% emgr - EMpirical GRamian Framework
%
%  project: emgr ( http://gramian.de )
%  version: 5.3-lgc ( 2018-01-01 )
%  authors: Christian Himpe ( 0000-0003-2194-6754 )
%  license: BSD 2-Clause License ( opensource.org/licenses/BSD-2-Clause )
%  summary: Empirical Gramians for (nonlinear) input-output systems.
%
% USAGE:
%
%  W = emgr_lgc(f,g,s,t,w,[pr],[nf],[ut],[us],[xs],[um],[xm],[dp])
%
% DESCRIPTION:
%
% Empirical gramian matrix and empirical covariance matrix computation
% for model reduction, decentralized control, nonlinearity quantification,
% sensitivity analysis, parameter identification, uncertainty quantification &
% combined state and parameter reduction of large-scale input-output systems.
% Data-driven analysis of input-output coherence and system-gramian-based
% nonlinear model order reduction. Compatible with OCTAVE and MATLAB.
%
% ARGUMENTS:
%
%   f {handle} vector field handle: x' = f(x,u,p,t)
%   g {handle} output function handle: y = g(x,u,p,t)
%   s {vector} system dimensions: [inputs,states,outputs]
%   t {vector} time discretization: [time-step,time-horizon]
%   w {string} single character encoding gramian type:
%    * 'c' empirical controllability gramian (Wc)
%    * 'o' empirical observability gramian (Wo)
%    * 'x' empirical cross gramian (Wx aka Wco or Xcg)
%    * 'y' empirical linear cross gramian (Wy)
%    * 's' empirical sensitivity gramian (Ws)
%    * 'i' empirical identifiability gramian (Wi)
%    * 'j' empirical joint gramian (Wj)
%  pr {matrix|0} parameters, each column is one set
%  nf {vector|0} option flags, twelve components, default zero:
%    * center: no(0), steady(1), last(2), mean(3), rms(4), midrange(5)
%    * input scales: single(0), linear(1), geometric(2), log(3), sparse(4)
%    * state scales: single(0), linear(1), geometric(2), log(3), sparse(4)
%    * input rotations: unit(0), single(1)
%    * state rotations: unit(0), single(1)
%    * normalizing (only: Wc, Wo, Wx, Wy): no(0), Jacobi(1), steady-state(2)
%    * cross gramian type (only: Wx, Wy, Wj): regular(0), non-symmetric(1)
%    * extra input (only: Wo, Wx, Ws, Wi, Wj): no(0), yes(1)
%    * parameter centering (only: Ws, Wi, Wj): no(0), linear(1), log(2)
%    * Schur-complement (only: Wi, Wj): detailed(0), approximate(1)
%    * cross gramian partition size (only: Wx, Wj): full(0), partitioned(<N)
%    * cross gramian partition index (only: Wx, Wj): partition(>0)
%  ut {handle|1} input function handle: u_t = ut(t), default: delta-impulse(1)
%  us {vector|0} steady-state input
%  xs {vector|0} steady-state and initial state x0
%  um {matrix|1} input scales
%  xm {matrix|1} initial-state scales
%  dp {handle|@mtimes} custom inner product handle: z = dp(x,y)
%
% RETURNS:
%
%  W {matrix} Gramian Matrix (for: Wc, Wo, Wx, Wy)
%  W  {cell}  [State-, Parameter-] Gramian (for: Ws, Wi, Wj)
%
% CITATION:
%
%  C. Himpe (2018). emgr - EMpirical GRamian Framework (Version 5.3)
%  [Software]. Available from http://gramian.de . doi: 10.5281/zenodo.1134429
%
% SEE ALSO:
%  gram
%
% KEYWORDS:
%
%  model reduction, system gramians, empirical gramians, cross gramian, MOR
%
% Further information: http://gramian.de

%% ARGUMENT SETUP

    global ODE; % Integrator Handle
    if(isa(ODE,'function_handle')==0), ODE = @ssp2; end;

    % Version Info
    if(strcmp(f,'version')), W = 5.3; return; end;

    % Default Arguments
    if( (nargin<6)  || isempty(pr) ), pr = 0.0; end;
    if( (nargin<7)  || isempty(nf) ), nf = 0.0; end;
    if( (nargin<8)  || isempty(ut) ), ut = 1.0; end;
    if( (nargin<9)  || isempty(us) ), us = 0.0; end;
    if( (nargin<10) || isempty(xs) ), xs = 0.0; end;
    if( (nargin<11) || isempty(um) ), um = 1.0; end;
    if( (nargin<12) || isempty(xm) ), xm = 1.0; end;
    if( (nargin<13) || isempty(dp) ), dp = @mtimes; end;

%% GENERAL SETUP

    % System Dimensions
    M = s(1);                   % Number of inputs
    N = s(2);                   % Number of states
    Q = s(3);                   % Number of outputs
    A = (numel(s)==4) * s(end); % Number of augmented parameter-states
    P = size(pr,1);             % Dimension of parameter
    K = size(pr,2);             % Number of parameter-sets
    h = t(1);                   % Width of time step
    L = floor(t(2)/h) + 1;      % Number of time steps plus initial value

    % Lazy Arguments
    if(isnumeric(g) && g==1)    % Assume unit output functional
        g = @id;
        Q = N;
    end;

    w = lower(w); % Ensure lower case gramian type

    nf = [nf(:)',zeros(1,12-numel(nf))]; % Ensure flag vector length

    if(isnumeric(ut) && numel(ut)==1) % Built-in input functions
        if(ut==Inf) % Linear Chirp Input
            mh = ones(M,1);
            sh = (1.0/L - 0.1/h) / L;
            ut = @(t) 0.5 + mh*0.5*cos(2.0*pi*(((0.1/h)+0.5*sh*t).*t));
        else        % Delta Impulse Input
            mh = ones(M,1)./h;
            ut = @(t) mh * (t<=h);
        end;
    end;

    if(numel(us)==1), us = us * ones(M,1); end;
    if(numel(xs)==1), xs = xs * ones(N,1); end;
    if(numel(um)==1), um = um * ones(M,1); end;
    if(numel(xm)==1), xm = xm * ones(N-(w=='y')*(N-M),1); end;

    if(size(um,2)==1), um = scales(um,nf(2),nf(4)); end;
    if(size(xm,2)==1), xm = scales(xm,nf(3),nf(5)); end;

    C = size(um,2); % Number of input scales sets
    D = size(xm,2); % Number of state scales sets

%% GRAMIAN SETUP

    if( (w=='c' || w=='o' || w=='x' || w=='y') && nf(6) && A==0) % Normalizing
        TX = ones(N,1);
        switch(nf(6))

            case 1 % Jacobi preconditioner
                NF = nf;
                NF(6) = 0;
                DP = @(x,y) sum(x.*y',2); % Diagonal-only pseudo-kernel
                WT = emgr_lgc(f,g,s,t,w,pr,NF,ut,us,xs,um,xm,DP);
                TX = sqrt(WT);

            case 2 % Steady-state and steady-state input
                TX(xs~=0) = xs(xs~=0);
        end
        tx = 1.0./TX;
        F = f; f = @(x,u,p,t) tx.*F(TX.*x,u,p,t);
        G = g; g = @(x,u,p,t)     G(TX.*x,u,p,t);
        xs = tx.*xs;
    end

    if( (w=='o' || w=='x' || w=='s') && nf(8) ) % Extra input
        up = @(t) us + ut(t);
    else
        up = @(t) us;
    end;

%% GRAMIAN COMPUTATION

    switch(w) % Empirical gramian types

        % Common layout:
        %  Setup: Initialize variable for empirical gramian matrix
        %  Loop nesting: for-each {parameter, input|state(parameter) scale}
        %  Loop bodies: perturb, simulate, center, normalize, store|accumulate
        %  Post-processing: normalize, (symmetrize), (decompose)
        %  Parameter-space gramians call state-space gramians

        case 'c' % Controllability gramian
            W = 0; % Reserve gramian variable
            for k = 1:K
                for c = 1:C
                    for m = find(um(:,c))' % parfor
                        em = sparse(m,1,um(m,c),M,1);
                        uu = @(t) us + ut(t) .* em;
                        x = ODE(f,@id,t,xs,uu,pr(:,k));
                        x = bsxfun(@minus,x,avg(x,nf(1),xs));
                        x = x * (1.0/um(m,c));
                        W = W + dp(x,x');
                    end;
                end;
            end;
            W = W * (h/(C*K));
            W = 0.5 * (W + W');

        case 'o' % Observability gramian
            W = 0; % Reserve gramian variable
            o = zeros(Q*L,N+A); % Pre-allocate observability matrix
            for k = 1:K
                for d = 1:D
                    for n = find(xm(:,d))' % parfor
                        xx = xs + sparse(n,1,xm(n,d),N+A,1);
                        if(A==0), pp = pr(:,k); else, pp = xx(N+1:end); end;
                        y = ODE(f,g,t,xx(1:N),up,pp);
                        y = bsxfun(@minus,y,avg(y,nf(1),g(xs(1:N),us,pp,0)));
                        y = y * (1.0/xm(n,d));
                        o(:,n) = y(:);
                    end;
                    W = W + dp(o',o);
                end;
            end;
            W = W * (h/(D*K));
            W = 0.5 * (W + W');

        case 'x' % Cross gramian
            assert(M==Q || nf(7),'emgr: non-square system!');

            i0 = 1;
            i1 = N+A;

            if(nf(11)>0) % Partitioned cross gramian

                i0 = i0 + (nf(12) - 1) * nf(11);
                i1 = min(i0 + nf(11) - 1,N);

                if(i0>N)
                    i0 = i0 - ( ceil( N / nf(11) ) * nf(11) - N);
                    i1 = min(i0 + nf(11) - 1,N+A);
                end;

                if(i0>i1 || i0<0), W = 0; return; end;
            end;

            W = 0; % Reserve gramian (partition) variable
            o = zeros(L,i1-i0+1,Q); % Pre-allocate observability 3-tensor

            for k = 1:K
                for d = 1:D
                    for n = find(xm(i0:i1,d))' % parfor
                        xx = xs + sparse(i0-1+n,1,xm(n,d),N+A,1);
                        if(A==0), pp = pr(:,k); else, pp = xx(N+1:end); end;
                        y = ODE(f,g,t,xx(1:N),up,pp);
                        y = bsxfun(@minus,y,avg(y,nf(1),g(xs(1:N),us,pp,0)));
                        y = y * (1.0/xm(n,d));
                        o(:,n,:) = y';
                    end;
                    if(nf(7)) % Non-symmetric cross gramian: cache average
                        o(:,:,1) = sum(o,3);
                    end;
                    for c = 1:C % parfor
                        for m = find(um(:,c))'
                            em = sparse(m,1,um(m,c),M,1);
                            uu = @(t) us + ut(t) .* em;
                            if(A==0), pp = pr(:,k); else, pp = xs(N+1:end); end;
                            x = ODE(f,@id,t,xs(1:N),uu,pp);
                            x = bsxfun(@minus,x,avg(x,nf(1),xs(1:N)));
                            x = x * (1.0/um(m,c));
                            if(nf(7)) % Non-symmetric cross gramian
                                W = W + dp(x,o(:,:,1));
                            else      % Regular cross gramian
                                W = W + dp(x,o(:,:,m));
                            end;
                        end;
                    end;
                end;
            end;
            W = W * (h/(C*D*K));

        case 'y' % Linear cross gramian
            assert(M==Q || nf(7),'emgr: non-square system!');

            W = 0; % Reserve gramian variable
            o = zeros(L,N,M); % Pre-allocate adjoint 3-tensor

            for k = 1:K
                for c = 1:C
                    for m = find(xm(:,c))' % parfor
                        em = sparse(m,1,xm(m,c),Q,1);
                        uu = @(t) us + ut(t) .* em;
                        z = ODE(g,@id,t,xs,uu,pr(:,k));
                        z = bsxfun(@minus,z,avg(z,nf(1),xs));
                        z = z * (1.0/xm(m,c));
                        o(:,:,m) = z';
                    end;
                    if(nf(7)) % Non-symmetric cross gramian: cache average
                        o(:,:,1) = sum(o,3);
                    end;
                    for m = find(um(:,c))' % parfor
                        em = sparse(m,1,um(m,c),M,1);
                        uu = @(t) us + ut(t) .* em;
                        x = ODE(f,@id,t,xs,uu,pr(:,k));
                        x = bsxfun(@minus,x,avg(x,nf(1),xs));
                        x = x * (1.0/um(m,c));
                        if(nf(7)) % Non-symmetric cross gramian
                            W = W + dp(x,o(:,:,1));
                        else      % Regular cross gramian
                            W = W + dp(x,o(:,:,m));
                        end;
                    end;
                end;
            end;
            W = W * (h/(C*K));

        case 's' % Sensitivity gramian
            [pr,pm] = pscales(pr,nf(9),size(um,2));
            W{1} = emgr_lgc(f,g,[M,N,Q],t,'c',pr,nf,ut,us,xs,um,xm,dp);
            W{2} = zeros(P,1); % Sensitivity gramian diagonal
            DP = @(x,y) sum(sum(x.*y')); % Trace pseudo-kernel
            UT = @(t) 1.0;
            for k = 1:P
                ek = sparse(k,1,1.0,P,1);
                F = @(x,u,p,t) f(x,up(t),pr + u * ek,t);
                W{2}(k) = emgr_lgc(F,g,[1,1,Q],t,'c',0,nf,UT,0,xs,pm(k,:),xm,DP);
            end;

        case 'i' % Identifiability gramian
            [pr,pm] = pscales(pr,nf(9),size(xm,2));
            V = emgr_lgc(f,g,[M,N,Q,P],t,'o',0,nf,ut,us,[xs;pr],um,[xm;pm],dp);
            W{1} = V(1:N,1:N); % Observability gramian
            WM = V(1:N,N+1:N+P);
            if(nf(10))         % Identifiability gramian
                W{2} = V(N+1:N+P,N+1:N+P);
            else
                W{2} = V(N+1:N+P,N+1:N+P) - (WM' * ainv(W{1}) * WM);
            end;

        case 'j' % Joint gramian
            [pr,pm] = pscales(pr,nf(9),size(xm,2));
            V = emgr_lgc(f,g,[M,N,Q,P],t,'x',0,nf,ut,us,[xs;pr],um,[xm;pm],dp);
            if(nf(11)) % Partitioned joint gramian
                W = V;
                return;
            end;
            W{1} = V(1:N,1:N); % Cross gramian
            WM = V(1:N,N+1:N+P);
            if(nf(10))         % Cross-identifiability gramian
                W{2} = -0.5 * (WM' * WM);
            else
                W{2} = -0.5 * (WM' * ainv(W{1} + W{1}') * WM);
            end;

        otherwise
            error('emgr: unknown gramian type!');
    end;
end

%% LOCALFUNCTION: scales
%  summary: Input and initial state perturbation scales
function sm = scales(s,d,c)

    switch(d)

        case 1 % Linear
            sc = [0.25,0.50,0.75,1.0];

        case 2 % Geometric
            sc = [0.125,0.25,0.5,1.0];

        case 3 % Logarithmic
            sc = [0.001,0.01,0.1,1.0];

        case 4 % Sparse
            sc = [0.01,0.5,0.99,1.0];

        otherwise % One
            sc = 1;
    end;

    if(c==0), sc = [-sc,sc]; end;

    sm = s * sc;
end

%% LOCALFUNCTION: pscales
%  summary: Parameter perturbation scales
function [pr,pm] = pscales(p,d,c)

    assert(size(p,2)>=2,'emgr: min + max parameter required!');

    pmin = min(p,[],2);
    pmax = max(p,[],2);

    switch(d) % Parameter centering

        case 1 % Linear
            pr = 0.5 * (pmax + pmin);
            pm = (pmax - pmin) * linspace(0,1.0,c);
            pm = bsxfun(@plus,pm,pmin - pr);

        case 2 % Logarithmic
            lmin = log(pmin);
            lmax = log(pmax);
            pr = real(exp(0.5 * (lmax + lmin)));
            pm = (lmax - lmin) * linspace(0,1.0,c);
            pm = bsxfun(@plus,pm,lmin);
            pm = bsxfun(@plus,real(exp(pm)),pmin - pr);

        otherwise % None
            pr = pmin;
            pm = (pmax - pmin) * linspace(0,1.0,c);
    end;
end

%% LOCALFUNCTION: id
%  summary: output identity function
function x = id(x,u,p,t)

end

%% LOCALFUNCTION: avg
%  summary: State and output trajectory centering
function mn = avg(x,d,c)

    switch(d)

        case 1 % Steady state / output
            mn = c;

        case 2 % Final state / output
            mn = x(:,end);

        case 3 % Mean state / output
            mn = mean(x,2);

        case 4 % Root-mean-square state / output
            mn = sqrt(mean(x.*x,2));

        case 5 % Midrange state / output
            mn = 0.5*(max(x,[],2)-min(x,[],2));

        otherwise % None
            mn = zeros(size(x,1),1);
    end;
end

%% LOCALFUNCTION: ainv
%  summary: Quadratic complexity approximate inverse matrix
function x = ainv(m)

    d = diag(m);
    d(d~=0) = 1.0./d(d~=0);
    x = bsxfun(@times,m,-d);
    x = bsxfun(@times,x,d');
    x(1:numel(d)+1:end) = d;
end

%% LOCALFUNCTION: ssp2
%  summary: Low-Storage Strong Stability Preserving Runge-Kutta SSPx2
function y = ssp2(f,g,t,x0,u,p)

    global STAGES;

    if(isscalar(STAGES)==0), STAGES = 3; end;

    h = t(1);
    K = floor(t(2)/h) + 1;

    y = g(x0,u(0),p,0);
    y(end,K) = 0; % Pre-allocate trajectory

    hs = h / (STAGES-1);
    xk1 = x0;
    xk2 = x0;

    for k = 2:K
        tk = (k - 1.5) * h;
        uk = u(tk);
        for s = 1:(STAGES-1)
            xk1 = xk1 + hs * f(xk1,uk,p,tk);
            tk = tk + hs;
        end;
        xk1 = (xk2 + (STAGES-1) * xk1 + h * f(xk1,uk,p,tk)) ./ STAGES;
        xk2 = xk1;
        y(:,k) = g(xk1,uk,p,tk);
    end;
end
