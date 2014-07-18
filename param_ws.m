function param_ws(o)
% param_ws (parameter reduction)
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
 X = ones(N,1);

 rand('seed',1009);
 A = rand(N,N); A(1:N+1:end) = -0.55*N;
 B = rand(N,J);
 C = rand(O,N);
 P = 0.01*rand(N,1);

 LIN = @(x,u,p) A*x+B*u+p;
 OUT = @(x,u,p) C*x;

%%%%%%%% Parameter Reduction %%%%%%%%%

% FULL
 tic;
 Y = ab2(LIN,OUT,T,X,U,P);
 FULL = toc

% OFFLINE
 tic;
 WS = emgr(LIN,OUT,[J N O],T,'s',P,0,1,0,X);
 [PP D QQ] = svd(WS{2}); PP = PP(1:R,:); QQ = QQ(1:R,:)';
 p = QQ*PP*P;
 OFFLINE = toc

% ONLINE
 tic;
 y = ab2(LIN,OUT,T,X,U,p);
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
 print -dsvg param_ws.svg;

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
