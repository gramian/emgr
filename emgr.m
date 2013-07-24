function W = emgr(f,g,q,p,t,w,nf,ut,us,xs,um,xm,yd)
% emgr - Empirical Gramian Framework ( Version: 1.2 )
% by Christian Himpe, 2013 ( http://gramian.de )
% released under BSD 2-Clause License ( opensource.org/licenses/BSD-2-Clause )
%
% SYNOPSIS:
%    W = emgr(f,g,q,p,t,w,[nf],[ut],[us],[xs],[um],[xm],[yd]);
%
% ABOUT:
%    emgr - Empirical Gramian Framemwork, computating empirical gramians
%    for model reduction and system identification.
%    Compatible with OCTAVE and MATLAB.
%
% INPUTS:
%   (func handle)  f - system function handle, signature: xdot = f(x,u,p)
%   (func handle)  g - output function handle, signature:    y = g(x,u,p)
%        (vector)  q - system dimensions [inputs,states,outputs]
%        (vector)  p - parameters, 0 if none
%        (vector)  t - time [start,step,stop]
%          (char)  w - gramian type:
%            * 'C' or 'c' : empirical controllability gramian (WC)
%            * 'O' or 'o' : empirical observability gramian (WO)
%            * 'X' or 'x' : empirical cross gramian (WX, also WCO)
%            * 'S' or 's' : empirical sensitivity gramian (WS)
%            * 'I' or 'i' : empirical identifiability gramian (WI)
%            * 'J' or 'j' : empirical joint gramian (WJ)
%        (vector) [nf = 0] - options, 10 components
%            + residual steady(0), mean(1), median(2), last(3), pod(4)
%            + unit-normal(0), pod(1) directions
%            + linear(0), log(1), geometric(2), single(3) input scale spacing
%            + linear(0), log(1), geometric(2), single(3) init-state scale spacing
%            + unit(0), [factorial(1)], dyadic(2), single(3) input rotations
%            + unit(0), [factorial(1)], dyadic(2), single(3) init-state rotations
%            + single(0), double(1) run
%            + disable(0), enable(1)
%                * robust parameters (WC, WS only)
%                * data-driven pod (WO, WI only)
%                * enforce symmetry (WX, WJ only)
%            + disable(0), enable(1) data-driven gramians
%            + solver: euler(0), adams-bashforth(1), leapfrog(2)
%  (matrix,vector,scalar) [ut = 1] - input; default: delta impulse
%         (vector,scalar) [us = 0] - steady-state input
%         (vector,scalar) [xs = 0] - steady-state
%  (matrix,vector,scalar) [um = 1] - input scales
%  (matrix,vector,scalar) [xm = 1] - init-state scales
%           (cell,matrix) [yd = 0] - observed data
%
% OUTPUT:
%            (matrix)  W - Gramian matrix (WC, WO, WX only)
%              (cell)  W - {State-,Parameter-} Gramian Matrices (WS, WI, WJ only)
%
% TODO:
%     factorial transformations
%
% For further information see http://gramian.de
%*

J = q(1);             % number of inputs
N = q(2);             % number of states
O = q(3);             % number of outputs
p = p(:);
P = numel(p);         % number of parameters
T = (t(3)-t(1))/t(2); % number of time steps
h = t(2);             % time step width
w = lower(w);

if (nargin<7) ||(isempty(nf)), nf = 0; end;
if (nargin<8) ||(isempty(ut)), ut = 1; end;
if (nargin<9) ||(isempty(us)), us = 0; end;
if (nargin<10)||(isempty(xs)), xs = 0; end;
if (nargin<11)||(isempty(um)), um = 1; end;
if (nargin<12)||(isempty(xm)), xm = 1; end;
if (nargin<13)||(isempty(yd)), yd = {sparse(N,T);sparse(O,T)}; end;

if (numel(nf)<10), nf(10) = 0; end;
if (numel(ut)==1), ut = ones(J,1)*ut; end;
if (numel(us)==1), us = ones(J,1)*us; end;
if (numel(xs)==1), xs = ones(N,1)*xs; end;
if (numel(um)==1), um = ones(J,1)*um; end;
if (numel(xm)==1), xm = ones(N,1)*xm; end;

if(w=='c' || w=='o' || w=='x')

    if(w=='c'&&nf(8)~=0)
        J = J+P;
        if(size(ut,1)==J-P), ut = [ut;ones(P,1)]; end;
        if(size(us,1)==J-P), us = [us;P*ones(P,1)]; end;
        if(size(um,1)==J-P), um = [um;ones(P,1)]; end;
        F = f; f = @(x,u,p) F(x,u(1:J-P),u(J-P+1:J));
        G = g; g = @(x,u,p) G(x,u(1:J-P),u(J-P+1:J));
    end

    if(size(ut,2)==1), ut = [ut,sparse(J,T-1)]; k = (1.0/h); else k = 1.0; end;
    if(size(us,2)==1), us = us*ones(1,T); end;
    if(size(um,2)==1), um = scales(um,nf(3),nf(5)); end;
    if(size(xm,2)==1), xm = scales(xm,nf(4),nf(6)); end;
    C = size(um,2);
    D = size(xm,2);

    if(nf(2)==1)&&(w~='o'), dx=svd(ut,'econ');                    else dx=0; end;
    if(nf(2)==1)&&(w~='c'), dy=svd(odex(f,N,h,T,xs,us,p),'econ'); else dy=0; end;
    if(nf(1)==0), X = xs; Y = g(xs,zeros(J,1),p); else X = 0; Y = 0; end;

    W = zeros(N,N);
    o = zeros(O,T,N);

    if(nf(7)==1) % double run
        nf(7) = 0;
        A = emgr(f,g,q,p,t,w,nf,ut,us,xs,um,xm,yd);
        A = sqrt(diag(A));
        B = diag(1.0./A);
        A = diag(A);

        F = f;
        G = g;
        f = @(x,u,p) A*F(B*x,u,p);
        g = @(x,u,p)   G(B*x,u,p);
    end

    if(nf(9)==1) % data driven
        if(size(yd,1)==1 && w=='o'), yd = {[];yd{:}}; end;
        C = size(yd,1); um = ones(J,C);
        D = size(yd,2); xm = ones(N,D);
    end
end

switch(w)

    case 'c' % controllability gramian
        for c=1:C
            for j=1:J % parfor
                uu = us + bsxfun(@times,ut,dirs(j,J,dx)*(um(j,c)*k));
                if(nf(9)==0), x=odex(f,h,T,xs,uu,p,nf(10)); else x=yd{1,c}; end;
                x = bsxfun(@minus,x,res(nf(1),x,X))*(1.0/um(j,c));
                W = W + x*x';
            end
        end
        W = W*(h/C);

    case 'o' % observability gramian
        for d=1:D
            for n=1:N % parfor
                xx = xs + dirs(n,N,dy)*xm(n,d);
                if(nf(9)==0), y=odey(f,g,h,T,xx,us,O,p,nf(10)); else y=yd{2,d}; end;
                o(:,:,n) = bsxfun(@minus,y,res(nf(1),y,Y))*(1.0/xm(n,d));
            end
            for n=1:N
                for m=1:N
                    W(n,m) = W(n,m) + sum(sum(o(:,:,n).*o(:,:,m)));
                end
            end
        end
        W = W*(h/D);

    case 'x' % cross gramian
        if(J~=O), error('ERROR: non-square system!'); end;
        for d=1:D
            for n=1:N % parfor
                xx = xs + dirs(n,N,dy)*xm(n,d);
                if(nf(9)==0), y=odey(f,g,h,T,xx,us,O,p,nf(10)); else y=yd{2,d}; end;
                o(:,:,n) = bsxfun(@minus,y,res(nf(1),y,Y))*(1.0/xm(n,d));
            end
            for c=1:C
                for j=1:J % parfor
                    uu = us + bsxfun(@times,ut,dirs(j,J,dx)*(um(j,c)*k));
                    if(nf(9)==0), x=odex(f,h,T,xs,uu,p,nf(10)); else x=yd{1,c}; end;
                    x = bsxfun(@minus,x,res(nf(1),x,X))*(1.0/um(j,c));
                    W = W + x*permute(o(j,:,:),[2 3 1]);
                end
            end
        end
        W = W*(h/(C*D));

    case 's' % sensitivity gramian
        W = cell(2,1);
        W{1} = emgr(f,g,[J N O],zeros(P,1),t,'c',nf,ut,us,xs,um,xm);
        W{2} = eye(P);
        F = @(x,u,p) f(x,zeros(J,1),p*u);
        G = @(x,u,p) g(x,zeros(J,1),p*u);
        for q=1:P
            V = emgr(F,G,[1 N O],(1:P==q),t,'c',nf,1,p(q),xs,1,xm);
            W{1} = W{1} + V;      % controllability gramian
            W{2}(q,q) = trace(V); % sensitivity gramian
        end

    case 'i' % identifiability gramian
        if(size(xm,1)==N), xm = [xm;ones(P,1)]; end;
        W = cell(2,1);
        F = @(x,u,p) [f(x(1:N),u,x(N+1:N+P));zeros(P,1)];
        G = @(x,u,p)  g(x(1:N),u,x(N+1:N+P));
        V = emgr(F,G,[J N+P O],0,t,'o',nf,ut,us,[xs;p],um,xm);
        W{1} = V(1:N,1:N);         % observability gramian
        W{2} = V(N+1:N+P,N+1:N+P); % identifiability gramian

    case 'j' % joint gramian
        if(size(xm,1)==N), xm = [xm;ones(P,1)]; end;
        if(size(um,1)==J), um = [um;ones(P,1)]; end;
        W = cell(2,1);
        F = @(x,u,p) [f(x(1:N),u(1:J),x(N+1:N+P));u(J+1:J+P)];
        G = @(x,u,p) [g(x(1:N),u(1:J),x(N+1:N+P));x(N+1:N+P)];
        V = emgr(F,G,[J+P N+P O+P],0,t,'x',nf,[ut;ones(P,1)],[us;zeros(P,1)],[xs;p],um,xm);
        W{1} = V(1:N,1:N);         % cross gramian
        W{2} = V(N+1:N+P,N+1:N+P); % cross identifiability gramian

    otherwise
        error('ERROR: unknown gramian type!');
end

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
    end

    switch(e)
        case 0 % unit
            s = [-s,s];
        case 1 % factorial
            s = s*(2*(dec2bin(0:2^q-1)-'0')'-1)./sqrt(2^q);
        case 2 % dyadic
            s = s*s';
        case 3 % single
            %s = s;
    end
end

%%%%%%%%
function d = dirs(n,N,e)

    switch(e)
        case 0    % unit-normal
            d = (1:N==n)';
        otherwise % POD
            d = e;
    end
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
    end
end

%%%%%%%%
function x = odex(f,h,T,Z,u,p,q)

z = Z;
x = zeros(numel(z),T);

    switch(q)
        case 0 % Eulers Method
            for t=1:T
                z = z + h*f(z,u(:,t),p);
                x(:,t) = z;
            end
        case 1 % Adams-Bashforth Method
            m = 0.5*h*f(z,u(:,1),p);
            z = z + h*f(z + m,u(:,1),p);
            x(:,1) = z;

            for t=2:T
                k = 0.5*h*f(z,u(:,t),p);
                z = z + 3.0*k - m;
                x(:,t) = z;
                m = k;
            end
        case 2 % Leapfrog Method
            N = size(z,1);
            m = N/2;
            n = m + 1;
            l = f(z,u(:,1),p);
            k = l(n:N);

            for t=1:T
                z(1:m) = z(1:m) + h*z(n:N) + 0.5*h*h*k;
                l = f(z,u(:,t),p);
                z(n:N) = z(n:N) + 0.5*h*(k+l(n:N));
                x(:,t) = z;
                k = l(n:N);
            end
    end
end

%%%%%%%%
function y = odey(f,g,h,T,Z,u,o,p,q)

z = Z;
y = zeros(o,T);

    switch(q)
        case 0 % Eulers Method
            for t=1:T
                z = z + h*f(z,u(:,t),p);
                y(:,t) = g(z,u(:,t),p);
            end
        case 1 % Adams-Bashforth Method
            m = 0.5*h*f(z,u(:,1),p);
            z = z + h*f(z + m,u(:,1),p);
            y(:,1) = g(z,u(:,1),p);

            for t=2:T
                k = 0.5*h*f(z,u(:,t),p);
                z = z + 3.0*k - m;
                y(:,t) = g(z,u(:,t),p);
                m = k;
            end
        case 2 % Leapfrog Method
            N = size(z,1);
            m = N/2;
            n = m + 1;
            l = f(z,u(:,1),p);
            k = l(n:N);

            for t=1:T
                z(1:m) = z(1:m) + h*z(n:N) + 0.5*h*h*k;
                l = f(z,u(:,t),p);
                z(n:N) = z(n:N) + 0.5*h*(k+l(n:N));
                y(:,t) = g(z,u(:,t),p);
                k = l(n:N);
            end
    end
end

