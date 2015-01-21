function nbody(o)
% nbody (n-body reduction)
% by Christian Himpe, 2013-2015 ( http://gramian.de )
% released under BSD 2-Clause License ( opensource.org/licenses/BSD-2-Clause )
%*

if(exist('emgr')~=2)
    disp('emgr framework is required. Download at http://gramian.de/emgr.m');
    return;
end

%% Setup

N = 5;
T = [0.0,0.01,2.0];
q = [0,4*N,2*N];
R = 4;
p = ones(N,1);

X = [1.449;  0.0;    0.400; -0.345; -1.125;...
     0.448; -1.125; -0.448;  0.400;  0.345;...
     0.0;   -0.922; -1.335;  0.810; -0.919;...
    -0.349;  0.919; -0.349;  1.335;  0.810]; %5-body eight-figure
F = @(x,u,p) [x((2*N)+1:end);acc(x(1:2*N),u,p)];
G = @(x,u,p)  x(1:2*N);

%% Main

Y = irk3(F,G,T,X,zeros(1,200),p); % Full Order

tic;
WO = emgr(F,G,q,T,'o',p);
WOP = WO(1:(2*N),1:(2*N));
WOV = WO((2*N)+1:end,(2*N)+1:end);
[PP DD QQ] = svd(WOP); PP = PP(:,1:2*R); QQ = PP'; %diag(DD)'
[TT DD VV] = svd(WOV); TT = TT(:,1:2*R); VV = TT'; %diag(DD)'
OFFLINE = toc


f = @(x,u,p) [x((2*R)+1:end);QQ*acc(PP*x(1:2*R),u,p)];
g = @(x,u,p) PP*x(1:2*R);
x = [QQ*X(1:2*N);QQ*X((2*N)+1:end)];
y = irk3(f,g,T,x,zeros(1,200),p);

%{
f = @(x,u,p) [x((2*R)+1:end);VV*acc(TT*x(1:2*R),u,p)];
g = @(x,u,p) TT*x(1:2*R);
x = [VV*X(1:2*N);VV*X((2*N)+1:end)];
y = irk3(f,g,T,x,zeros(1,200),p);
%}

%%%%%%%% Plot %%%%%%%%

if(nargin==0), return; end
figure();
hold on;
cmap = hsv(N+1);
d = 1;
for c=1:N
    plot(y(d,:),y(d+1,:),'--','Color',cmap(c,:));
    plot(y(d,end),y(d+1,end),'*','Color',cmap(c,:));
    d = d + 2;
end
hold off;
ylim([-0.6 0.6]);
pbaspect([2,1,1]);
set(gcf,'InvertHardcopy','off')
if(o==1), print('-dpng',[mfilename(),'.png']); end;

%% ======== Acceleration ========

function y = acc(x,u,p)

    N = numel(x)/2;
    A = reshape(x,[2 N]);
    B = zeros(size(A));
    Z = zeros(1,N);
    y = zeros(2,N);

    for I=1:N
        B = bsxfun(@minus,A,A(:,I));
        Z = p'./(sqrt(0.0001+(sum(B.^2))).^3);
        B = bsxfun(@times,B,Z);
        y(:,I) = sum(B,2);
    end

    y = y(:);

%% ======== Leapfrog ========

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
