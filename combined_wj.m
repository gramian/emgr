function combined_wj(o)
% combined nonlinear reduction
% by Christian Himpe, 2013-2014 ( http://gramian.de )
% released under BSD 2-Clause License ( opensource.org/licenses/BSD-2-Clause )
%*

if(exist('emgr')~=2) disp('emgr framework is required. Download at http://gramian.de/emgr.m'); return; end

%%%%%%%% Setup %%%%%%%%

 J = 8;
 O = J;
 N = 64;
 R = O;
 T = [0 0.01 1];
 L = (T(3)-T(1))/T(2);
 U = [ones(J,1) zeros(J,L-1)];
 X = zeros(N,1);

 rand('seed',1009);
 A = rand(N,N); A(1:N+1:end) = -0.55*N; A = 0.5*(A+A');
 B = rand(N,J);
 C = B';
 P = rand(N,1);

 NON = @(x,u,p) A*tanh(p.*x) + B*u;
 OUT = @(x,u,p) C*x;

%%%%%%%% Reduction %%%%%%%%

% FULL
 tic;
 Y = ab2(NON,OUT,T,X,U,P);
 FULL = toc

% OFFLINE
 tic;
 WJ = emgr(NON,OUT,[J N O],T,'j',P);
 [UU D VV] = svd(WJ{1}); UU = UU(:,1:R);   VV = UU';
 [PP D QQ] = svd(WJ{2}); PP = PP(1:R*R,:); QQ = PP';
 x = VV*X;
 p = PP*P;
 non = @(x,u,p) VV*NON(UU*x,u,QQ*p);
 out = @(x,u,p) OUT(UU*x,u,QQ*p);
 OFFLINE = toc

% ONLINE
 tic;
 y = ab2(non,out,T,x,U,p);
 ONLINE = toc

%%%%%%%% Output %%%%%%%%

% TERMINAL
 ERROR = norm(norm(Y - y)./norm(Y))
 RELER = abs(Y - y)./abs(Y);

% PLOT 
 l = (1:-0.01:0)'; cmap = [l,l,ones(101,1)]; cmax = max(max(RELER));
 figure('PaperSize',[2.4,6.4],'PaperPosition',[0,0,6.4,2.4]);
 imagesc(RELER); caxis([0 cmax]); cbr = colorbar; colormap(cmap); 
 set(gca,'YTick',1:N,'xtick',[]); set(cbr,'YTick',[0 cmax],'YTickLabel',{'0',sprintf('%0.1e',cmax)});
 if(nargin>0), print -dsvg combined_wj.svg; end

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
