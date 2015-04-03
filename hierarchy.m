function hierarchy(o)
% hierarchy (hierarchical network reduction)
% by Christian Himpe, 2013-2015 ( http://gramian.de )
% released under BSD 2-Clause License ( opensource.org/licenses/BSD-2-Clause )
%*

if(exist('emgr')~=2)
    disp('emgr framework is required. Download at http://gramian.de/emgr.m');
    return;
end

%% Setup

D = 3;		%Tree depth
M = 4;		%Children per node

J = 1;
O = M^D;
N = (M^(D+1)-1)/(M-1);
T = [0,1,100];
h = T(2);
L = (T(3)-T(1))/T(2);
U = exp(-0.0005*(1:L).^2);
X = zeros(N,1);

rand('seed',1009);
A = trasm(D,M);
B = sparse(N,1); B(1,1) = D;
C = [sparse(O,N-O),speye(O)];

LIN = @(x,u,p) A*x + B*u;
OUT = @(x,u,p) C*x;

ADJ = @(x,u,p) A'*x + C'*u;
AOU = @(x,u,p) B'*x;

norm1 = @(y) T(2)*sum(abs(y(:)));
norm2 = @(y) sqrt(T(2)*dot(y(:),y(:)));
norm8 = @(y) max(y(:));

%% Main

Y = C*rk1(LIN,h,L,X,U,0); % Full Order

global CUSTOM_ODE;
CUSTOM_ODE = @rk1;

tic;
WC = emgr(LIN,OUT,[J,N,O],T,'c',0,[0,0,0,0,0,0,0,0,0,0,0,-1]);
WO = emgr(ADJ,AOU,[O,N,J],T,'c',0,[0,0,0,0,0,0,0,0,0,0,0,-1]);
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
    y = c*rk1(lin,h,L,x,U,0); % Reduced Order
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

%% ======== Tree Assembler ========

function A = trasm(d,c)

 a = (c^(d+1)-1)/(c-1);
 A = -speye(a);

 for I=0:((a-1)/c)-1
  b = 1+c*I;
  %A(1+I,b+1:b+c) = rand(1,c)+1;
  A(b+1:b+c,1+I) = rand(c,1)+1;
 end

%% ======== Balancer ========

function [X Y Z] = balance(WC,WO)

 [L D l] = svd(WC); LC = L*diag(sqrt(diag(D)));
 [L D l] = svd(WO); LO = L*diag(sqrt(diag(D)));
 [U Y V] = svd(LO'*LC);
 X = ( LO*U*diag(1.0./sqrt(diag(Y))) )';
 Z =   LC*V*diag(1.0./sqrt(diag(Y)));

%% ======== Integrator ========

function y = rk1(f,h,T,x,u,p)

    y(:,1) = x;
    y(end,T) = 0;

    for t=1:T
        x = h*f(x,u(:,t),p);
        y(:,t) = x;
    end;
