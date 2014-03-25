function benchmark_non(o)
% benchmark (nonlinear rc ladder)
% by Christian Himpe, 2013-2014 ( http://gramian.de )
% released under BSD 2-Clause License ( opensource.org/licenses/BSD-2-Clause )
%*

if(exist('emgr')~=2) disp('emgr framework is required. Download at http://gramian.de/emgr.m'); return; end

%%%%%%%% Setup %%%%%%%%

 J = 1;
 N = 64;
 O = 1;
 R = 8;
 T = [0 0.01 1.0];
 L = (T(3)-T(1))/T(2);
 U = ones(1,L);
 X = zeros(N,1);

 g = @(x) exp(x)+x-1.0;
 
 A1 = spdiags(ones(N-1,1),-1,N,N)-speye(N);
 A2 = spdiags([ones(N-1,1);0],0,N,N)-spdiags(ones(N,1),1,N,N);
 
 NON = @(x,u,p) g(A1*x)-g(A2*x) + [u;sparse(N-1,1)];
 OUT = @(x,u,p) x(1);

%%%%%%%% Reduction %%%%%%%%%

% FULL
 tic;
 Y = ab2(NON,OUT,T,X,U,0);
 FULL = toc

% OFFLINE
 tic;
 WX = emgr(NON,OUT,[J N O],T,'x');
 [UU D VV] = svd(WX); UU = UU(:,1:R); VV = VV(:,1:R)';
 x = VV*X;
 non = @(x,u,p) VV*NON(UU*x,u,p);
 out = @(x,u,p) OUT(UU*x,u,p);
 OFFLINE = toc

% ONLINE
 tic;
 y = ab2(non,out,T,x,U,0);
 ONLINE = toc

%%%%%%%% Output %%%%%%%%

% TERMINAL
 ERROR = norm(norm(Y - y)./norm(Y))
 RELER = abs(Y - y)./abs(Y);

% PLOT 
 l = (1:-0.01:0)'; cmap = [l,l,ones(101,1)]; cmax = max(max(RELER));
 figure('PaperSize',[2.4,6.4],'PaperPosition',[0,0,6.4,2.4]);
 imagesc([RELER;RELER;RELER]); caxis([0 cmax]); cbr = colorbar; colormap(cmap); 
 set(gca,'YTick',1:N,'Yticklabel',{'','1',''},'xtick',[]); set(cbr,'YTick',[0 cmax],'YTickLabel',{'0',sprintf('%0.1e',cmax)});
 if(nargin>0), print -dsvg benchmark_non.svg; end

%%%%%%%% Integrator %%%%%%%%

function y = ab2(f,g,t,x,u,p)

 h = t(2);
 T = t(3)/h;
 m = f(x,u(:,1),p);
 x = x + (0.5*h)*(m + f(x + h*m,u(:,1),p));
 y(:,1) = g(x,u(:,1),p);

 for t=2:T
     k = (0.5*h)*f(x,u(:,t),p);
     x = x + 3.0*k - m;
     y(:,t) = g(x,u(:,1),p);
     m = k;
 end;

