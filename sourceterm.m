function sourceterm(o)
% sourceterm reduction
% by Christian Himpe, 2013 ( http://gramian.de )
% released under BSD 2-Clause License ( http://gramian.de/#license )
%*

if(exist('emgr')~=2) disp('emgr framework is required. Download at http://gramian.de/emgr.m'); return; end

%%%%%%%% Setup %%%%%%%%

 J = 10;
 O = J;
 N = 100;
 R = 10;
 T = [0 0.01 1];
 L = (T(3)-T(1))/T(2);
 U = [ones(J,1) zeros(J,L-1)];
 X =  ones(N,1);

 A = rand(N,N);
 A(1:N+1:end) = -0.5*N;
 B = rand(N,J);
 C = rand(O,N);
 P = 0.1*rand(N,1);

 LIN = @(x,u,p) A*x+B*u+p;
 OUT = @(x,u,p) C*x;

%%%%%%%% Parameter Reduction %%%%%%%%%

% FULL
 tic; Y = rk2(LIN,OUT,[J N O],T,X,U,P); ORIGINAL = toc

% OFFLINE
 tic;
 WS = emgr(LIN,OUT,[J N O],P,T,'s',1,1,0,X);
 [PP D QQ] = svd(WS{2}); PP = PP(1:R,:); QQ = QQ(1:R,:)';
 p = QQ*PP*P;
 OFFLINE = toc

% ONLINE
 tic; y = rk2(LIN,OUT,[J N O],T,X,U,p); ONLINE = toc

%%%%%%%% Output %%%%%%%%

% TERMINAL
 ERROR = norm(norm(Y - y)./norm(Y))
 RELER = abs(Y - y)./abs(Y);

% PLOT
 if(nargin<1 || o==0 ), return; end
 l = (1:-0.01:0)'; cmap = [l,l,ones(101,1)];
 figure('PaperSize',[2.4,6.4],'PaperPosition',[0,0,6.4,2.4]);
 imagesc(RELER); caxis([0 max(max(RELER))]); colorbar; colormap(cmap); set(gca,'YTick',1:N);
 if(o==2 && exist('OCTAVE_VERSION')), print -dsvg source.svg; end

%%%%%%%% Integrator %%%%%%%%

function y = rk2(f,g,q,t,x,u,p)

 T = (t(3)-t(1))/t(2);
 y = zeros(q(3),T);
 h = t(2);

 for A=1:T
  x = x + h*f(x + 0.5*h*f(x,u(:,A),p),u(:,A),p); %Improved Eulers Method
  y(:,A) = g(x,u(:,A),p);
 end
