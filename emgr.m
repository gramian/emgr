function W = emgr(f,g,q,t,w,pr,nf,ut,us,xs,um,xm,yd)
% emgr - Empirical Gramian Framework ( Version: 1.5 )
% by Christian Himpe, 2013 ( http://gramian.de )
% released under BSD 2-Clause License ( opensource.org/licenses/BSD-2-Clause )
%
% SYNOPSIS:
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
%            * 's' : empirical sensitivity gramian (WS)
%            * 'i' : empirical identifiability gramian (WI)
%            * 'j' : empirical joint gramian (WJ)
%        (vector) [pr = 0] - parameters
%        (vector) [nf = 0] - options, 10 components:
%            + residual steady(0), mean(1), median(2), last(3), pod(4)
%            + unit-normal(0), pod(1) directions
%            + linear(0), log(1), geometric(2), single(3) input scale spacing
%            + linear(0), log(1), geometric(2), single(3) state scale spacing
%            + unit(0), [factorial(1)], dyadic(2), single(3) input rotations
%            + unit(0), [factorial(1)], dyadic(2), single(3) state rotations
%            + single(0), double(1) run
%            + disable(0), enable(1)
%                * robust parameters (WC, WS only)
%                * data-driven pod (WO, WI only)
%                * enforce symmetry (WX, WJ only)
%            + disable(0), enable(1) data-driven gramians
%            + solver: Euler(0), Adams-Bashforth(1), Leapfrog(2)
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
% KEYWORDS:
%    model reduction, empirical gramian, emgr
%
% TODO:
%     factorial transformations
%
% For further information see <http://gramian.de>
%*

w = lower(w);

J = q(1);             % number of inputs
N = q(2);             % number of states
O = q(3);             % number of outputs

M = N;
if(numel(q)==4), M = q(4); end;

T = (t(3)-t(1))/t(2); % number of time steps
h = t(2);             % time step width

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

if (numel(nf)<10), nf(10) = 0; end;
if (numel(ut)==1), ut(1:J,1) = ut; end;
if (numel(us)==1), us(1:J,1) = us; end;
if (numel(xs)==1), xs(1:N,1) = xs; end;
if (numel(um)==1), um(1:J,1) = um; end;
if (numel(xm)==1), xm(1:N,1) = xm; end;

if(w=='c' || w=='o' || w=='x')

    if(nf(7)==1) % double run
        nf(7) = 0;
        S = sqrt(diag(emgr(f,g,q,t,w,p,nf,ut,us,xs,um,xm,yd)));
        A = spdiags(S,0,N,N);
        B = spdiags(1.0./S,0,N,N);

        F = f;
        G = g;
        f = @(x,u,p) A*F(B*x,u,p);
        g = @(x,u,p)   G(B*x,u,p);
    end;

    if(w=='c'&&nf(8)~=0) % robust parameters
        J = J+P;
        if(size(us,1)==J-P), us = [us;p]; end;
        if(size(ut,1)==J-P), ut = [ut;ones(P,1)]; end;
        if(size(um,1)==J-P), um = [um;ones(P,1)]; end;
        F = f; f = @(x,u,p) F(x,u(1:J-P),u(J-P+1:J));
        G = g; g = @(x,u,p) G(x,u(1:J-P),u(J-P+1:J));
    end;

    if(size(ut,2)==1), ut(:,2:T) = 0; k = (1.0/h); else k = 1.0; end;
    if(size(us,2)==1), us = repmat(us,[1 T]); end;
    if(size(um,2)==1), um = scales(um,nf(3),nf(5)); end;
    if(size(xm,2)==1), xm = scales(xm,nf(4),nf(6)); end;
    C = size(um,2);
    D = size(xm,2);

    X = xs(1:M); 
    Y = g(xs(1:M),us,p); 

    dx = 0;
    dy = 0;
    if(nf(2)==1),
         dx = svd(ut,'econ');
         dy = svd(ode(f,h,T,xs,us,p,cf(10)),'econ');
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

    case 'c' % controllability gramian
        for c=1:C
            for j=1:J % parfor
                uu = us + bsxfun(@times,ut,dirs(j,J,dx)*(um(j,c)*k));
                if(nf(9)~=0), x = yd{1,c}; else
                    x = ode(f,h,T,xs,uu,p,nf(10));
                end;
                x = bsxfun(@minus,x,res(nf(1),x,X))*(1.0/um(j,c));
                W = W + x*x';
            end;
        end;
        W = W*(h/C);

    case 'o' % observability gramian
        for d=1:D
            for n=1:N % parfor
                xx = xs + dirs(n,N,dy)*xm(n,d);
                pp = p;
                if(nf(9)~=0), y = yd{2,d}; else
                    if(M~=N), pp = xx(M+1:N); xx = xx(1:M); end;
                    x = ode(f,h,T,xx,us,pp,nf(10));
                    for s=1:T, y(:,s) = g(x(:,s),us(:,1),pp); end;
                end;
                o(:,:,n) = bsxfun(@minus,y,res(nf(1),y,Y))*(1.0/xm(n,d));
            end;
            for n=1:N
                for m=1:N
                    W(n,m) = W(n,m) + sum(sum(o(:,:,n).*o(:,:,m)));
                end;
            end;
        end;
        W = W*(h/D);

    case 'x' % cross gramian
        if(J~=O), error('ERROR: non-square system!'); end;
        for d=1:D
            for n=1:N % parfor
                xx = xs + dirs(n,N,dy)*xm(n,d);
                pp = p;
                if(nf(9)~=0), y = yd{2,d}; else
                    if(M~=N), pp = xx(M+1:N); xx = xx(1:M); end;
                    x = ode(f,h,T,xx,us,pp,nf(10));
                    for s=1:T, y(:,s) = g(x(:,s),us(:,1),pp); end;
                end;
                o(:,:,n) = bsxfun(@minus,y,res(nf(1),y,Y))*(1.0/xm(n,d));
            end;
            for c=1:C
                for j=1:J % parfor
                    uu = us + bsxfun(@times,ut,dirs(j,J,dx)*(um(j,c)*k));
                    xx = xs;
                    pp = p;
                    if(nf(9)~=0), x = yd{1,c}; else
                        if(M~=N), pp = xs(M+1:N); xx = xs(1:M); end;
                        x = ode(f,h,T,xx,uu,pp,nf(10));
                    end;
                    x = bsxfun(@minus,x,res(nf(1),x,X))*(1.0/um(j,c));
                    W = W(1:M,:) + x*permute(o(j,:,:),[2 3 1]);
                end;
            end;
        end;
        W = W*(h/(C*D));

    case 's' % sensitivity gramian
        W = cell(2,1);
        W{1} = emgr(f,g,[J N O],t,'c',sparse(P,1),nf,ut,us,xs,um,xm);
        W{2} = eye(P); % speye
        F = @(x,u,p) f(x,us,p*u);
        G = @(x,u,p) g(x,us,p*u);
        for q=1:P
            V = emgr(F,G,[1 N O],t,'c',(1:P==q),nf,1,p(q),xs,1,xm);
            W{1} = W{1} + V;      % controllability gramian
            W{2}(q,q) = trace(V); % sensitivity gramian
        end;

    case 'i' % identifiability gramian
        if(size(xm,1)==N), xm = [xm;ones(P,1)]; end;
        W = cell(2,1);
        V = emgr(f,g,[J N+P O N],t,'o',p,nf,ut,us,[xs;p],um,xm);
        W{1} = V(1:N,1:N);         % observability gramian
        W{2} = V(N+1:N+P,N+1:N+P); % identifiability gramian

    case 'j' % joint gramian
        if(size(xm,1)==N), xm = [xm;ones(P,1)]; end;
        W = cell(2,1);
        V = emgr(f,g,[J N+P O N],t,'x',p,nf,ut,us,[xs;p],um,xm);
        W{1} = V(1:N,1:N);                       % cross gramian
        U = spdiags((1.0/diag(W{1}))',0,N,N);
        % U = U - U*(W{1}-diag(diag(W{1})))*U;
        W{2} = V(1:N,N+1:N+P)'*U*V(1:N,N+1:N+P); % cross identifiability gramian

    otherwise
        error('ERROR: unknown gramian type!');
end;

if(w=='c' || w=='o' || (w=='x'&&nf(8)==1)), W = 0.5*(W+W'); end;

end

%%%%%%%%
function s = scales(s,d,e)

    switch(d)
        case 0 % linear
            s = s*[0.25,0.50,0.75,1.0];
        case 1 % logarithmic
            s = s*[0.001,0.01,0.1,1.0];
        case 2 % geometric
            s = s*[0.125,0.25,0.5,1.0];
        case 3 % single
            %s = s;
    end;

    switch(e)
        case 0 % unit
            s = [-s,s];
        case 1 % factorial
            s = s*(2*(dec2bin(0:2^q-1)-'0')'-1)./sqrt(2^q);
        case 2 % dyadic
            s = s*s';
        case 3 % single
            %s = s;
    end;
end

%%%%%%%%
function d = dirs(n,N,e)

    switch(e)
        case 0    % unit-normal
            d = (1:N==n)';
        otherwise % POD
            d = e;
    end;
end

%%%%%%%%
function y = res(v,d,e)

    switch(v)
        case 0 % steady
            y = e;
        case 1 % average
            y = mean(d,2);
        case 2 % median
            y = median(d,2);
        case 3 % last
            y = d(:,end);
        case 4 % POD
            y = svd(d,'econ');
    end;
end

%%%%%%%%
function x = ode(f,h,T,z,u,p,O)

x(numel(z),T) = 0;

switch(O)

    case 0 % Eulers Method
        for t=1:T
            z = z + h*f(z,u(:,t),p);
            x(:,t) = z;
        end;

    case 1 % Adams-Bashforth Method
        m = 0.5*h*f(z,u(:,1),p);
        z = z + h*f(z + m,u(:,1),p);
        x(:,1) = z;

        for t=2:T
            k = 0.5*h*f(z,u(:,t),p);
            z = z + 3.0*k - m;
            x(:,t) = z;
            m = k;
        end;

    case 2 % Leapfrog Method
        N = numel(z);
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
end;

end

