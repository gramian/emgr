function benchmark_ilp(o)
% benchmark (inverse lyapunov procedure)
% by Christian Himpe, 2013-2014 ( http://gramian.de )
% released under BSD 2-Clause License ( opensource.org/licenses/BSD-2-Clause )
%*

if(exist('emgr')~=2) disp('emgr framework is required. Download at http://gramian.de/emgr.m'); return; end
if(exist('ilp')~=2)  disp('ilp generator is required. Download at http://gramian.de/ilp.m'); return; end

%%%%%%%% Setup %%%%%%%%

 J = 8;
 N = 64;
 O = J;
 R = 3*O;
 [A,B,C] = ilp(J,N,O,0,1009);
 T = [0.0,0.01,1.0];
 L = (T(3)-T(1))/T(2);
 U = [ones(J,1),zeros(J,L-1)];
 X =  zeros(N,1);

 LIN = @(x,u,p) A*x + B*u;
 OUT = @(x,u,p) C*x;

%%%%%%%% Reduction %%%%%%%%%

% FULL
 tic;
 Y = ab2(LIN,OUT,T,X,U,0);
 FULL = toc

% OFFLINE
 tic;
 WC = emgr(LIN,OUT,[J,N,O],T,'c');
 WO = emgr(LIN,OUT,[J,N,O],T,'o');
 [UU D VV] = balance(WC,WO,R);
 a = UU*A*VV;
 b = UU*B;
 c = C*VV;
 x = UU*X;
 lin = @(x,u,p) a*x + b*u;
 out = @(x,u,p) c*x;
 OFFLINE = toc

% ONLINE
 tic;
 y = ab2(lin,out,T,x,U,0);
 ONLINE = toc

%%%%%%%% Output %%%%%%%%

% TERMINAL
 norm2 = @(y) sqrt(T(2)*sum(sum(y.*y)));
 ERROR = norm2(Y - y)./norm2(Y)
 RELER = abs(Y - y)./abs(Y);

% PLOT
 if(nargin==0), return; end
 l = (1:-0.01:0)'; cmap = [l,l,ones(101,1)]; cmax = max(max(RELER));
 figure('PaperSize',[2.4,6.4],'PaperPosition',[0,0,6.4,2.4]);
 imagesc(RELER); caxis([0 cmax]); cbr = colorbar; colormap(cmap);
 set(gca,'YTick',1:N,'xtick',[]); set(cbr,'YTick',[0 cmax],'YTickLabel',{'0',sprintf('%0.1e',cmax)});
 print -dsvg benchmark_ilp.svg;

%%%%%%%% Balancer %%%%%%%%

function [X Y Z] = balance(WC,WO,R)

 L = chol(WC+eye(size(WC,1)))-eye(size(WC,1));
 [U Y V] = svd(L*WO*L');
 X = diag(sqrt(diag(Y(1:R,1:R)))) * V(:,1:R)' / L';
 Z = L'*U(:,1:R)*diag(1./sqrt(diag(Y(1:R,1:R))));

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
