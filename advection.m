function advection(o)
% finite difference discretized pde transport equation reduction
% by Christian Himpe, 2013-2015 ( http://gramian.de )
% released under BSD 2-Clause License ( opensource.org/licenses/BSD-2-Clause )
%*

if(exist('emgr')~=2)
    disp('emgr framework is required. Download at http://gramian.de/emgr.m');
    return;
end

%% Setup

N = 100;
J = 0;
O = N;
R = N/10;
T = [0.0,0.01,1.0];
L = (T(3)-T(1))/T(2);
U = zeros(1,L);
H = 0.1;
X = 0:H:9.9; X = exp(-(X-1).^2)';

P = 5;
A = spdiags((1.0/H)*[ones(N,1) -ones(N,1)],[-1,0],N,N);

LIN = @(x,u,p) p*A*x;
ADJ = @(x,u,p) p*A'*x + u;
OUT = @(x,u,p) x;

%% Main

Y = irk3(LIN,OUT,T,X,U,P); % Full Order

tic;
WO = emgr(ADJ,1,[N,N,O],T,'c',P);
[UU D VV] = svd(WO); UU = UU(:,1:R); VV = VV(:,1:R)';
a = VV*A*UU;
c = UU;
x = VV*X;
lin = @(x,u,p) p*a*x;
out = @(x,u,p) c*x;
OFFLINE = toc

y = irk3(lin,out,T,x,U,P);

%% Output

if(nargin==0), return; end
cmax = max(y(:));
figure();
imagesc(sparse(y)); caxis([0 cmax]);
colormap([(1:-0.01:0)',(1:-0.01:0)',ones(101,1)]);
set(gca,'YTick',0,'xtick',[]); ylabel('X');
daspect([1,1,1]);
if(o==1), print('-dsvg',[mfilename(),'.svg']); end;

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
