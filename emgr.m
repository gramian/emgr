function W = emgr(f,g,q,t,w,pr,nf,ut,us,xs,um,xm,yd)
% emgr - Empirical Gramian Framework ( Version: 1.6 )
% by Christian Himpe 2013, 2014 ( http://gramian.de )
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
%            * 'y' : fast linear cross gramian (WY)
%            * 's' : empirical sensitivity gramian (WS)
%            * 'i' : empirical identifiability gramian (WI)
%            * 'j' : empirical joint gramian (WJ)
%        (vector) [pr = 0] - parameters
%        (vector) [nf = 0] - options, 11 components:
%            + residual steady(0), mean(1), median(2), last(3), pod(4)
%            + unit-normal(0), pod(1), factorial(2) directions
%            + linear(0), log(1), geometric(2), single(3) input scale spacing
%            + linear(0), log(1), geometric(2), single(3) state scale spacing
%            + unit(0), inverse(1), dyadic(2), single(3) input rotations
%            + unit(0), inverse(1), dyadic(2), single(3) state rotations
%            + single(0), double(1), scaled(2) run
%            + disable(0), enable(1) 
%                * robust parameters (WC, WS only)
%                * data-driven pod   (WO, WI only)
%                * enforce symmetry  (WX, WJ only)
%            + disable(0), enable(1) data-driven gramians
%            + solver: Euler(0), Adams-Bashforth(1), Leapfrog(2), Custom(-1)
%            + disable(0), enable(1) parameter scaling (WS, WI, WJ only)
%  (matrix,vector,scalar) [ut = 1] - input; default: delta impulse
%         (vector,scalar) [us = 0] - steady-state input
%         (vector,scalar) [xs = 0] - steady-state
%  (matrix,vector,scalar) [um = 1] - input scales
%  (matrix,vector,scalar) [xm = 1] - init-state scales
%           (cell,matrix) [yd = 0] - observed data
%
% OUTPUT:
%            (matrix)  W - Gramian matrix (WC, WO, WX only)
%              (cell)  W - {State-,Parameter-} Gramian Matrix (WS, WI, WJ only)
%
% SEE ALSO:
%    gram
%
% KEYWORDS:
%    model reduction, empirical gramian, emgr
%
%
% For further information see <http://gramian.de>
%*

w = lower(w);

J = q(1);             % number of inputs
N = q(2);             % number of states
O = q(3);             % number of outputs

M = N;
if(numel(q)==4), M = q(4); end;

h = t(2);             % time step width
T = (t(3)-t(1))/h;    % number of time steps

if (isnumeric(g) && g==1), g = @(x,u,p) x; O = N; end;

if (nargin<6) ||(isempty(pr)), pr = 0; end;
if (nargin<7) ||(isempty(nf)), nf = 0; end;
if (nargin<8) ||(isempty(ut)), ut = 1; end;
if (nargin<9) ||(isempty(us)), us = 0; end;
if (nargin<10)||(isempty(xs)), xs = 0; end;
if (nargin<11)||(isempty(um)), um = 1; end;
if (nargin<12)||(isempty(xm)), xm = 1; end;
if (nargin<13)||(isempty(yd)), yd = 0; end;

P = numel(pr);        % number of parameters
p = pr(:);

if (numel(nf)<11), nf(11)    = 0;  end;
if (numel(ut)==1), ut(1:J,1) = ut; end;
if (numel(us)==1), us(1:J,1) = us; end;
if (numel(xs)==1), xs(1:N,1) = xs; end;
if (numel(um)==1), um(1:J,1) = um; end;
if (numel(xm)==1), xm(1:N,1) = xm; end;

if(w=='c' || w=='o' || w=='x' || w=='y')

    switch(nf(7))
        case 1, % double run
            nf(7) = 0;
            TX = sqrt(diag(emgr(f,g,q,t,w,p,nf,ut,us,xs,um,xm,yd)));
            tx = 1.0./TX;

            F = f; f = @(x,u,p) TX.*F(tx.*x,u,p);
            G = g; g = @(x,u,p)     G(tx.*x,u,p);

        case 2, % scaled run
            TX = xs; tx = TX; tx(tx~=0) = 1.0./tx(tx~=0);
            TU = us; tu = TU; tu(tu~=0) = 1.0./tu(tu~=0);

            F = f; f = @(x,u,p) TX.*F(tx.*x,tu.*u,p);
            G = g; g = @(x,u,p)     G(tx.*x,tu.*u,p); 
    end;

    if(w=='c'&&nf(8)==1) % robust parameters
        J = J+P;
        if(size(us,1)==J-P), us = [us;p]; end;
        if(size(ut,1)==J-P), ut = [ut;ones(P,1)]; end;
        if(size(um,1)==J-P), um = [um;ones(P,1)]; end;
        F = f; f = @(x,u,p) F(x,u(1:J-P),u(J-P+1:J));
        G = g; g = @(x,u,p) G(x,u(1:J-P),u(J-P+1:J));
    end;

    X = xs(1:M);
    Y = g(xs(1:M),us(:,1),p);
    if(w=='y'), Y = X; end;

    if(size(ut,2)==1), ut(:,2:T) = 0; k = (1.0/h); else k = 1.0; end;
    if(size(us,2)==1), us = repmat(us,[1 T]); end;
    if(size(um,2)==1), um = scales(um,nf(3),nf(5)); end;
    if(size(xm,2)==1), xm = scales(xm,nf(4),nf(6)); end;
    C = size(um,2);
    D = size(xm,2);

    ud = 0;
    xd = 0;

    switch(nf(2)) % directions
        case 0, % unit-normal
            ud = eye(J);
            xd = eye(N);
        case 1, % pod
            [ud,E,V] = svd(ut,'econ');
            [xd,E,V] = svd(ode(f,h,T,xs,us,p,cf(10)),'econ');
        case 2, % factorial
            ud = (dec2bin(1:2^J-1)-'0')'*(2.0/sqrt(2^J)); % JJ = size(ud,2);
            xd = (dec2bin(1:2^N-1)-'0')'*(2.0/sqrt(2^N)); % NN = size(xd,2);
    end;

    if(nf(9)==1) % data driven
        if(size(yd,1)==1 && w=='o'), yd = {[];yd{:}}; end;
        C = size(yd,1); um = ones(J,C);
        D = size(yd,2); xm = ones(N,D);
    end;

    y = zeros(O,T);
    o = zeros(O,T,N);
    W = zeros(N,N);
end;

switch(w)

    case 'c', % controllability gramian
        for c=1:C
            for j=1:J % parfor
                uu = us + bsxfun(@times,ut,ud(:,j)*(um(j,c)*k));
                if(nf(9)==1), x = yd{1,c}; else
                    x = ode(f,h,T,xs,uu,p,nf(10));
                end;
                x = bsxfun(@minus,x,res(nf(1),x,X))*(1.0/um(j,c));
                W = W + x*x';
            end;
        end;
        W = W*(h/C);

    case 'o', % observability gramian
        for d=1:D
            for n=1:N % parfor
                xx = xs + xd(:,n)*xm(n,d);
                pp = p;
                if(nf(9)==1), y = yd{2,d}; else
                    if(M~=N), pp = xx(M+1:N); xx = xx(1:M); end;
                    x = ode(f,h,T,xx,us,pp,nf(10));
                    for s=1:T, y(:,s) = g(x(:,s),us(:,1),pp); end;
                end;
                o(:,:,n) = bsxfun(@minus,y,res(nf(1),y,Y))*(1.0/xm(n,d));
            end;
            for n=1:N
                for m=1:n
                    W(n,m) = W(n,m) + sum(sum(o(:,:,n).*o(:,:,m)));
                    W(m,n) = W(n,m);
                end;
            end;
        end;
        W = W*(h/D);

    case 'x', % cross gramian
        if(J~=O), error('ERROR: (emgr) non-square system!'); end;
        for d=1:D
            for n=1:N % parfor
                xx = xs + xd(:,n)*xm(n,d);
                pp = p;
                if(nf(9)==1), y = yd{2,d}; else
                    if(M~=N), pp = xx(M+1:N); xx = xx(1:M); end;
                    x = ode(f,h,T,xx,us,pp,nf(10));
                    for s=1:T, y(:,s) = g(x(:,s),us(:,1),pp); end;
                end;
                o(:,:,n) = bsxfun(@minus,y,res(nf(1),y,Y))*(1.0/xm(n,d));
            end;
            for c=1:C
                for j=1:J % parfor
                    uu = us + bsxfun(@times,ut,ud(:,j)*(um(j,c)*k));
                    xx = xs;
                    pp = p;
                    if(nf(9)==1), x = yd{1,c}; else
                        if(M~=N), pp = xs(M+1:N); xx = xs(1:M); end;
                        x = ode(f,h,T,xx,uu,pp,nf(10));
                    end;
                    x = bsxfun(@minus,x,res(nf(1),x,X))*(1.0/um(j,c));
                    W = W(1:M,:) + x*permute(o(j,:,:),[2 3 1]);
                end;
            end;
        end;
        W = W*(h/(C*D));

    case 'y', % fast cross gramian
        if(J~=O), error('ERROR: (emgr) non-square system!'); end;
        for c=1:C
            for j=1:J % parfor
                uu = us + bsxfun(@times,ut,ud(:,j)*(um(j,c)*k));
                yy = us + bsxfun(@times,ut,ud(:,j)*(xm(j,c)*k));
                if(nf(9)==1), x = yd{1,c}; y = yd{2,c}; else
                    x = ode(f,h,T,xs,uu,p,nf(10));
                    y = ode(g,h,T,xs,yy,p,nf(10));
                end;
                x = bsxfun(@minus,x,res(nf(1),x,X))*(1.0/um(j,c));
                y = bsxfun(@minus,y,res(nf(1),y,Y))*(1.0/xm(j,c));
                W = W + x*y';
            end;
        end;
        W = W*(h/C);

    case 's', % sensitivity gramian
        W = cell(2,1);
        W{1} = emgr(f,g,[J N O],t,'c',sparse(P,1),nf,ut,us,xs,um,xm);
        W{2} = eye(P);
        F = @(x,u,p) f(x,us,p*u);
        G = @(x,u,p) g(x,us,p*u);
        for q=1:P
            up = 1; if(nf(11)~=0), up = p(q); end;
            V = emgr(F,G,[1 N O],t,'c',(1:P==q),nf,1,p(q),xs,up,xm);
            W{1} = W{1} + V;      % controllability gramian
            W{2}(q,q) = trace(V); % sensitivity gramian
        end;

    case 'i', % identifiability gramian
        xp = ones(P,1); if(nf(11)~=0), xp = p; end; 
        if(size(xm,1)==N), xm = [xm;xp]; end;
        W = cell(2,1);
        V = emgr(f,g,[J N+P O N],t,'o',p,nf,ut,us,[xs;p],um,xm);
        W{1} = V(1:N,1:N);         % observability gramian
        W{2} = V(N+1:N+P,N+1:N+P); % identifiability gramian

    case 'j', % joint gramian
        xp = ones(P,1); if(nf(11)~=0), xp = p; end;
        if(size(xm,1)==N), xm = [xm;xp]; end;
        W = cell(2,1);
        V = emgr(f,g,[J N+P O N],t,'x',p,nf,ut,us,[xs;p],um,xm);
        W{1} = V(1:N,1:N);                       % cross gramian
        U = spdiags((1.0/diag(W{1}))',0,N,N);
        % U = U - U*(W{1}-diag(diag(W{1})))*U;
        W{2} = V(1:N,N+1:N+P)'*U*V(1:N,N+1:N+P); % cross-identifiability gramian

    otherwise
        error('ERROR: (emgr) unknown gramian type!');
end;

if(w=='c' || w=='o' || ((w=='x'||w=='y')&&nf(8)==1)), W = 0.5*(W+W'); end;

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

%%%%%%%% RESIDUAL SELECTOR %%%%%%%%
function y = res(v,d,e)

    switch(v)
        case 0, % steady
            y = e;
        case 1, % average
            y = mean(d,2);
        case 2, % median
            y = median(d,2);
        case 3, % last
            y = d(:,end);
        case 4, % POD
            [U,E,V] = svd(d,'econ'); y = U(1,:);
    end;
end

%%%%%%%% INTEGRATOR %%%%%%%%
function x = ode(f,h,T,z,u,p,O)

N = numel(z);
x(N,T) = 0;

switch(O)

    case 0, % Eulers Method
        for t=1:T
            z = z + h*f(z,u(:,t),p);
            x(:,t) = z;
        end;

    case 1, % Adams-Bashforth Method
        m = 0.5*h*f(z,u(:,1),p);
        z = z + h*f(z + m,u(:,1),p);
        x(:,1) = z;

        for t=2:T
            k = 0.5*h*f(z,u(:,t),p);
            z = z + 3.0*k - m;
            x(:,t) = z;
            m = k;
        end;

    case 2, % Leapfrog Method
        n = N/2;
        l = f(z,u(:,1),p);
        k = l(n+1:N);

        for t=1:T
            z(1:n) = z(1:n) + h*z(n+1:N) + 0.5*h*h*k;
            l = f(z,u(:,t),p);
            z(n+1:N) = z(n+1:N) + 0.5*h*(k+l(n+1:N));
            x(:,t) = z;
            k = l(n+1:N);
        end;

    case -1, % Custom Solver
        x = CUSTOM_ODE(f,h,T,z,u,p,O);

end;

end

