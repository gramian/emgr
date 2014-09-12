function advection(o)
% finite difference discretized pde transport equation reduction
% by Christian Himpe, 2013-2014 ( http://gramian.de )
% released under BSD 2-Clause License ( opensource.org/licenses/BSD-2-Clause )
%*

if(exist('emgr')~=2) disp('emgr framework is required. Download at http://gramian.de/emgr.m'); return; end

%%%%%%%% Setup %%%%%%%%

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

%%%%%%%% Reduction %%%%%%%%%

% FULL
 tic;
 Y = ab2(LIN,OUT,T,X,U,P);
 FULL = toc

% OFFLINE
 tic;
 WO = emgr(ADJ,1,[N,N,O],T,'c',P);
 [UU D VV] = svd(WO); UU = UU(:,1:R); VV = VV(:,1:R)';
 a = VV*A*UU;
 c = UU;
 x = VV*X;
 lin = @(x,u,p) p*a*x;
 out = @(x,u,p) c*x;
 OFFLINE = toc

% ONLINE
 tic;
 y = ab2(lin,out,T,x,U,P);
 ONLINE = toc

%%%%%%%% Output %%%%%%%

% TERMINAL
 norm2 = @(y) sqrt(T(2)*sum(sum(y.*y)));
 ERROR = norm2(Y - y)./norm2(Y)

% PLOT
 if(nargin==0), return; end
 l = (1:-0.01:0)'; cmap = [l,l,ones(101,1)];
 figure('PaperSize',[2.4,2.6],'PaperPosition',[0,0,2.6,2.4]);
 imagesc(sparse(y)); colormap(cmap);
 set(gca,'YTick',0,'xtick',[]); ylabel('X');
 print -dpng advection.png;

%%%%%%%% Integrator %%%%%%%%

function y = ab2(f,g,t,x,u,p)

 h = t(2);
 T = t(3)/h;
 m = 0.5*h*f(x,u(:,1),p);
 x = x + h*f(x + m,u(:,1),p);
 y(:,1) = g(x,u(:,1),p);
 y(end,T) = 0;

 for t=2:T
     k = (0.5*h)*f(x,u(:,t),p);
     x = x + 3*k - m;
     y(:,t) = g(x,u(:,1),p);
     m = k;
 end;

