function benchmark_ilp(o)
% benchmark (inverse lyapunov procedure)
% by Christian Himpe, 2013-2015 ( http://gramian.de )
% released under BSD 2-Clause License ( opensource.org/licenses/BSD-2-Clause )
%*

if(exist('emgr')~=2)
    disp('emgr framework is required. Download at http://gramian.de/emgr.m');
    return;
end

if(exist('ilp')~=2)
    disp('ilp generator is required. Download at http://gramian.de/ilp.m');
    return;
end

%% Setup

J = 8;
N = 64;
O = J;
[A,B,C] = ilp(J,N,O,0,1009);
T = [0.0,0.01,1.0];
L = (T(3)-T(1))/T(2);
U = [ones(J,1),zeros(J,L-1)];
X =  zeros(N,1);

LIN = @(x,u,p) A*x + B*u;
OUT = @(x,u,p) C*x;

norm1 = @(y) T(2)*sum(abs(y(:)));
norm2 = @(y) sqrt(T(2)*dot(y(:),y(:)));
norm8 = @(y) max(y(:));

%% Main

Y = irk3(LIN,OUT,T,X,U,0); % Full Order

tic;
WC = emgr(LIN,OUT,[J,N,O],T,'c');
WO = emgr(LIN,OUT,[J,N,O],T,'o');
[VV D UU] = balance(WC,WO);
OFFLINE = toc

for I=1:N-1
    uu = UU(:,1:I);
    vv = VV(1:I,:);
    a = vv*A*uu;
    b = vv*B;
    c = C*uu;
    x = vv*X;
    lin = @(x,u,p) a*x + b*u;
    out = @(x,u,p) c*x;
    y = irk3(lin,out,T,x,U,0); % Reduced Order
    l1(I) = norm1(Y-y)/norm1(Y);
    l2(I) = norm2(Y-y)/norm2(Y);
    l8(I) = norm8(Y-y)/norm8(Y);
end;

%% Output

if(nargin==0), return; end
figure();
semilogy(1:N-1,l1,'r','linewidth',2); hold on;
semilogy(1:N-1,l2,'g','linewidth',2);
semilogy(1:N-1,l8,'b','linewidth',2); hold off;
xlim([1,N-1]);
ylim([10^floor(log10(min([l1(:);l2(:);l8(:)]))-1),1]);
pbaspect([2,1,1]);
legend('L1 Error ','L2 Error ','L8 Error ','location','northeast');
if(o==1), print('-dsvg',[mfilename(),'.svg']); end;

%% ======== Balancer ========

function [X Y Z] = balance(WC,WO)

 [L D l] = svd(WC); LC = L*diag(sqrt(diag(D)));
 [L D l] = svd(WO); LO = L*diag(sqrt(diag(D)));
 [U Y V] = svd(LO'*LC);
 X = ( LO*U*diag(1.0./sqrt(diag(Y))) )';
 Z =   LC*V*diag(1.0./sqrt(diag(Y)));

%% ======== Integrator ========

function y = irk3(f,g,t,x,u,p)

    h = t(2);
    T = round(t(3)/h);

    k1 = h*f(x,u(:,1),p);
    k2 = h*f(x + 0.5*k1,u(:,1),p);
    k3r = h*f(x + 0.75*k2,u(:,1),p);
    x = x + (2.0/9.0)*k1 + (1.0/3.0)*k2 + (4.0/9.0)*k3r; % Ralston RK3

    y(:,1) = g(x,u(:,1),p);
    y(end,T) = 0;

    for t=2:T
        l1 = h*f(x,u(:,t),p);
        l2 = h*f(x + 0.5*l1,u(:,t),p);
        x = x + (2.0/3.0)*l1 + (1.0/3.0)*k1 + (5.0/6.0)*(l2 - k2);
        y(:,t) = g(x,u(:,t),p);
        k1 = l1;
        k2 = l2;
    end;
