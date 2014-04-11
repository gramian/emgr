function W = emgr(f,g,q,t,w,pr,nf,ut,us,xs,um,xm,yd)
% emgr - Empirical Gramian Framework ( Version: 2.0 )
% by Christian Himpe 2013-2014 ( http://gramian.de )
% released under BSD 2-Clause License ( opensource.org/licenses/BSD-2-Clause )
%
% SYNTAX:
%    W = emgr(f,g,q,t,w,[pr],[nf],[ut],[us],[xs],[um],[xm],[yd]);
%
% ABOUT:
%    emgr - Empirical Gramian Framemwork, computating empirical gramians
%    for model reduction and system identification.
%    Compatible with OCTAVE and MATLAB.
%
% INPUTS:
%   (func handle)  f - system function handle; signature: xdot = f(x,u,p)
%   (func handle)  g - output function handle; signature:    y = g(x,u,p)
%        (vector)  q - system dimensions [inputs,states,outputs]
%        (vector)  t - time discretization [start,step,stop]
%          (char)  w - gramian type:
%            * 'c' : empirical controllability gramian (WC)
%            * 'o' : empirical observability gramian (WO)
%            * 'x' : empirical cross gramian (WX)
%            * 'y' : empirical approximate cross gramian (WY)
%            * 's' : empirical sensitivity gramian (WS)
%            * 'i' : empirical identifiability gramian (WI)
%            * 'j' : empirical joint gramian (WJ)
%        (vector) [pr = 0] - parameters
%        (vector) [nf = 0] - options, 12 components:
%            + residual steady(0), mean(1), median(2), last(3), pod(4)
%            + linear(0), log(1), geometric(2), single(3) input scale spacing
%            + linear(0), log(1), geometric(2), single(3) state scale spacing
%            + unit(0), inverse(1), dyadic(2), single(3) input rotations
%            + unit(0), inverse(1), dyadic(2), single(3) state rotations
%            + single(0), double(1), scaled(2) run
%            + default(0), data-driven gramians(1)
%            + default(0), active(1) passive(2) robust params; only: WC,WS,WY
%            + default(0), parameter scaling(1); only WS,WI,WJ
%            + default(0), schur complement(1) for param extract; only: WI
%            + default(0), enforce symmetry(1)
%            + Euler(0),Two-Step(1),Leapfrog(2),Ralston(3),Custom(-1) solver
%  (matrix,vector,scalar) [ut = 1] - input; default: delta impulse
%         (vector,scalar) [us = 0] - steady-state input
%         (vector,scalar) [xs = 0] - steady-state
%  (matrix,vector,scalar) [um = 1] - input scales
%  (matrix,vector,scalar) [xm = 1] - init-state scales
%           (cell,matrix) [yd = 0] - observed data
%
% OUTPUT:
%            (matrix)  W - Gramian Matrix (WC, WO, WX only)
%              (cell)  W - {State-,Parameter-} Gramian (WS, WI, WJ only)
%
% SEE ALSO:
%    gram
%
% KEYWORDS:
%    model reduction, empirical gramian, emgr
%
% Further information: <http://gramian.de>
%*
    w = lower(w);

    J = q(1); % number of inputs
    N = q(2); % number of states
    O = q(3); % number of outputs
    M = N; if(numel(q)==4), M = q(4); end;

    if (isnumeric(g) && g==1), g = @(x,u,p) x; O = N; end;

    h = t(2);                 % time step width
    T = round((t(3)-t(1))/h); % number of time steps

    if(nargin<6) ||(isempty(pr)), pr = 0.0; end;
    if(nargin<7) ||(isempty(nf)), nf = 0.0; end;
    if(nargin<8) ||(isempty(ut)), ut = 1.0/h; end;
    if(nargin<9) ||(isempty(us)), us = 0.0; end;
    if(nargin<10)||(isempty(xs)), xs = 0.0; end;
    if(nargin<11)||(isempty(um)), um = 1.0; end;
    if(nargin<12)||(isempty(xm)), xm = 1.0; end;
    if(nargin<13)||(isempty(yd)), yd = 0.0; end;

    P = numel(pr); % number of parameters
    p = pr(:);

    if (isa(ut,'function_handle')),
        uf = ut;
        ut = zeros(J,T);
        for l=1:T, ut(:,l) = uf(l*h); end;
    end;

    if(numel(nf)<12), nf(12)    = 0;  end;
    if(numel(ut)==1), ut(1:J,1) = ut; end;
    if(numel(us)==1), us(1:J,1) = us; end;
    if(numel(xs)==1), xs(1:N,1) = xs; end;
    if(numel(um)==1), um(1:J,1) = um; end;
    if(numel(xm)==1), xm(1:N,1) = xm; end;

    if(w=='c' || w=='o' || w=='x' || w=='y')

        switch(nf(6))
            case 1, % double run
                nf(6) = 0;
                TX = sqrt(diag(emgr(f,g,q,t,w,p,nf,ut,us,xs,um,xm,yd)));
                tx = 1.0./TX;
                F = f; f = @(x,u,p) TX.*F(tx.*x,u,p);
                G = g; g = @(x,u,p)     G(tx.*x,u,p);
            case 2, % scaled run
                TX = xs; tx = TX; tx(tx==0) = 1.0; tx = 1.0./tx;
                TU = us; tu = TU; tu(tu==0) = 1.0; tu = 1.0./tu;
                F = f; f = @(x,u,p) TX.*F(tx.*x,tu.*u,p);
                G = g; g = @(x,u,p)     G(tx.*x,tu.*u,p); 
        end;

        if(size(ut,2)==1), ut(:,2:T) = 0.0; end;

        switch(nf(8))
            case 1, % active parameters
                J = J+P;
                if(size(us,1)==J-P), us = [us;p]; end;
                if(size(ut,1)==J-P), ut = [ut;ones(P,T)]; end;
                if(size(um,1)==J-P), um = [um;ones(P,1)]; end;
                F = f; f = @(x,u,p) F(x,u(1:J-P),u(J-P+1:J));
                G = g; g = @(x,u,p) G(x,u(1:J-P),u(J-P+1:J));
            case 2, % passive parameters
                if(size(um,1)==J), um = [um;ones(P,1)]; end;
        end;

        X = xs(1:M);
        Y = g(xs(1:M),us(:,1),p); if(w=='y'), Y = X; end;

        if(size(us,2)==1), us = repmat(us,[1 T]); end;
        if(size(um,2)==1), um = scales(um,nf(2),nf(4)); end;
        if(size(xm,2)==1), xm = scales(xm,nf(3),nf(5)); end;
        C = size(um,2); % number of input scales
        D = size(xm,2); % number of state scales

        switch(nf(1)) % residuals
            case 0, % steady state
                res = @(d,e) e;
            case 1, % mean state
                res = @(d,e) mean(d,2);
            case 2, % median state
                res = @(d,e) median(d,2);
            case 3, % final state
                res = @(d,e) d(:,end);
            case 4, % principal direction
                res = @(d,e) prd(d);
        end;

        if(nf(7)==1) % data driven
            if(size(yd,1)==1 && w=='o'), yd = {[];yd{:}}; end;
            C = size(yd,1);
            D = size(yd,2);
            um = ones(J,C);
            xm = ones(N,D);
        end;

        if(w=='o'||w=='x')
            y = zeros(O,T);
            o = zeros(O,T,N);
        end;

        W = zeros(N,N); % preallocate gramian
    end;

    switch(w)

        case 'c', % controllability gramian
            for c=1:C
                for j=1:J % parfor
                    if(nf(7)==1), x = yd{1,c}; else
                        uu = us + bsxfun(@times,ut,(1:J==j)'*um(j,c));
                        pp = p; if(nf(8)==2), pp = p.*um(J+1:end,c); end;
                        x = ode(f,h,T,xs,uu,pp,nf(12));
                    end;
                    x = bsxfun(@minus,x,res(x,X))*(1.0/um(j,c));
                    W = W + x*x';
                end;
            end;
            W = W*(h/C);

        case 'o', % observability gramian
            for d=1:D
                for n=1:N % parfor
                    if(nf(7)==1), y = yd{2,d}; else
                        xx = xs + (1:N==n)'*xm(n,d);
                        pp = p; if(M<N), pp = xx(M+1:N); end;
                        x = ode(f,h,T,xx(1:M),us,pp,nf(12));
                        for s=1:T, y(:,s) = g(x(:,s),us(:,1),pp); end;
                    end;
                    o(:,:,n) = bsxfun(@minus,y,res(y,Y))*(1.0/xm(n,d));
                end;
                for n=1:N
                    for m=1:n
                        W(n,m) = W(n,m) + sum(dot(o(:,:,n),o(:,:,m)));
                        W(m,n) = W(n,m);
                    end;
                end;
            end;

            W = W*(h/D);

        case 'x', % cross gramian
            if(J~=O),error('ERROR! emgr: non-square system!');end;
            for d=1:D
                for n=1:N % parfor
                    if(nf(7)==1), y = yd{2,d}; else
                        xx = xs + (1:N==n)'*xm(n,d);
                        pp = p; if(M<N), pp = xx(M+1:N); end;
                        x = ode(f,h,T,xx(1:M),us,pp,nf(12));
                        for s=1:T, y(:,s) = g(x(:,s),us(:,1),pp); end;
                    end;
                    o(:,:,n) = bsxfun(@minus,y,res(y,Y))*(1.0/xm(n,d));
                end;
                for c=1:C
                    for j=1:J % parfor
                        if(nf(7)==1), x = yd{1,c}; else
                            uu = us + bsxfun(@times,ut,(1:J==j)'*um(j,c));
                            pp = p; if(M<N), pp = xs(M+1:N); end;
                            x = ode(f,h,T,xs(1:M),uu,pp,nf(12));
                        end;
                        x = bsxfun(@minus,x,res(x,X))*(1.0/um(j,c));
                        W = W(1:M,:) + x*permute(o(j,:,:),[2 3 1]);
                    end;
                end;
            end;
            W = W*(h/(C*D));

        case 'y', % approximate cross gramian
            if(J~=O && nf(8)==0),error('ERROR! emgr: non-square system!');end;
            for c=1:C
                for j=1:J % parfor
                    if(nf(7)==1), x = yd{1,c}; y = yd{2,c}; else
                        uu = us + bsxfun(@times,ut,(1:J==j)'*um(j,c));
                        yy = us + bsxfun(@times,ut,(1:J==j)'*xm(j,c));
                        pp = p; if(nf(8)==2), pp = p.*um(J+1:end,c); end;
                        x = ode(f,h,T,xs,uu,pp,nf(12));
                        y = ode(g,h,T,xs,yy,pp,nf(12));
                    end;
                    x = bsxfun(@minus,x,res(x,X))*(1.0/um(j,c));
                    y = bsxfun(@minus,y,res(y,Y))*(1.0/xm(j,c));
                    W = W + x*y';
                end;
            end;
            W = W*(h/C);

        case 's', % sensitivity gramian
            W = cell(2,1);
            W{1} = emgr(f,g,[J N O],t,'c',zeros(P,1),nf,ut,us,xs,um,xm);
            W{2} = eye(P);
            F = @(x,u,p) f(x,us,p*u);
            G = @(x,u,p) g(x,us,p*u);
            for q=1:P
                if(nf(9)==0), pm = 1.0; else pm = p(q); end;
                if(size(um,1)==J+P), pm = um(q,:); end;
                V = emgr(F,G,[1 N O],t,'c',(1:P==q),nf,1.0,p(q),xs,pm,xm);
                W{1} = W{1} + V;      % approximate controllability gramian
                W{2}(q,q) = trace(V); % sensitivity gramian
            end;

        case 'i', % identifiability gramian
            if(nf(9)==0), pm = ones(P,1); else pm = p; end; 
            if(size(xm,1)==N), xm = [xm;pm]; end;
            W = cell(2,1);
            V = emgr(f,g,[J N+P O N],t,'o',p,nf,ut,us,[xs;p],um,xm);
            W{1} = V(1:N,1:N);         % observability gramian
            W{2} = V(N+1:N+P,N+1:N+P); % approximate identifiability gramian
            if(nf(10)==1),
                W{2} = W{2} - V(N+1:N+P,1:N)*fastinv(W{1})*V(1:N,N+1:N+P);
            end;

        case 'j', % joint gramian
            if(nf(9)==0), pm = ones(P,1); else pm = p; end;
            if(size(xm,1)==N), xm = [xm;pm]; end;
            W = cell(2,1);
            V = emgr(f,g,[J N+P O N],t,'x',p,nf,ut,us,[xs;p],um,xm);
            W{1} = V(1:N,1:N); % cross gramian
            %W{2} = zeros(P);   % cross-identifiability gramian
            W{2} = -0.5*V(1:N,N+1:N+P)'*fastinv(W{1}+W{1}')*V(1:N,N+1:N+P);

        otherwise,
            error('ERROR! emgr: unknown gramian type!');
    end;

    if(nf(11)==1 && (w=='c'||w=='o'||w=='x'||w=='y') ), W = 0.5*(W+W'); end;
end

%%%%%%%% SCALE SELECTOR %%%%%%%%
function s = scales(s,d,e)

    switch(d)
        case 0, % linear
            s = s*[0.25,0.50,0.75,1.0];
        case 1, % logarithmic
            s = s*[0.001,0.01,0.1,1.0];
        case 2, % geometric
            s = s*[0.125,0.25,0.5,1.0];
        case 3, % single
            %s = s;
    end;

    switch(e)
        case 0, % unit
            s = [-s,s];
        case 1, % anti
            s = [1.0./s,s];
        case 2, % dyadic
            s = s*s';
        case 3, % single
            %s = s;
    end;
end

%%%%%%%% PRINCIPAL DIRECTION %%%%%%%%
function x = prd(d)

    [U,D,V] = svd(d,'econ');
    x = U(:,1);
end

%%%%%%%% FAST INVERSION %%%%%%%%
function x = fastinv(m)

    d = diag(m);
    n = numel(d);
    x = speye(n);
    x(1:n+1:end) = 1.0./d(:);
    x = x - x*(m-diag(d))*x;
end

%%%%%%%% INTEGRATORS %%%%%%%%
function x = ode(f,h,T,z,u,p,O)

    N = numel(z);
    x(N,T) = 0.0;

    switch(O)

        case 0, % Eulers Method
            for t=1:T
                z = z + h*f(z,u(:,t),p);
                x(:,t) = z;
            end;

        case 3, % Ralstons Method
            for t=1:T
                k = h*f(z,u(:,t),p);
                z = z + 0.25*k + (0.75*h)*f(z + (2.0/3.0)*k,u(:,t),p);
                x(:,t) = z;
            end;

        case 1, % Adams-Bashforth Method
            m = (0.5*h)*f(z,u(:,1),p);
            z = z + h*f(z + m,u(:,1),p);
            x(:,1) = z;

            for t=2:T
                k = (0.5*h)*f(z,u(:,t),p);
                z = z + 3.0*k - m;
                x(:,t) = z;
                m = k;
            end;

        case 2, % Leapfrog Method
            n = N/2;
            l = f(z,u(:,1),p);
            k = l(n+1:N);

            for t=1:T
                z(1:n) = z(1:n) + h*z(n+1:N) + (0.5*h*h)*k;
                l = f(z,u(:,t),p);
                z(n+1:N) = z(n+1:N) + (0.5*h)*(k+l(n+1:N));
                x(:,t) = z;
                k = l(n+1:N);
            end;

        case -1, % Custom Solver
            % global CUSTOM_ODE; % Uncomment to use as handle
            x = CUSTOM_ODE(f,h,T,z,u,p);
    end;
end

