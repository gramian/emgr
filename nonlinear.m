function nonlinear(o)
% nonlinear reduction
% by Christian Himpe, 2013,2014 ( http://gramian.de )
% released under BSD 2-Clause License ( opensource.org/licenses/BSD-2-Clause )
%*

if(exist('emgr')~=2) disp('emgr framework is required. Download at http://gramian.de/emgr.m'); return; end

%%%%%%%% Setup %%%%%%%%

 J = 4;
 O = J;
 N = 16;
 R = O;
 T = [0 0.01 1];
 L = (T(3)-T(1))/T(2);
 U = [N*ones(J,1) zeros(J,L-1)];
 X =    ones(N,1);

 rand('seed',1009);
 A = rand(N,N); A(1:N+1:end) = -0.55*N; A = 0.5*(A+A');
 B = rand(N,J);
 C = B';
 P = A(:);

 NON = @(x,u,p) reshape(p,[N N])*sinc(10*x) + B*u;
 OUT = @(x,u,p) C*x;

%%%%%%%% Reduction %%%%%%%%

% FULL
 FULL = cputime;
 Y = rk2(NON,OUT,T,X,U,P);
 FULL = cputime - FULL

% OFFLINE
 OFFLINE = cputime;
 WJ = emgr(NON,OUT,[J N O],T,'j',P,0,1,0,X);
 [UU D VV] = svd(WJ{1}); UU = UU(:,1:R);   VV = VV(:,1:R)';
 [PP D QQ] = svd(WJ{2}); PP = PP(1:R*R,:); QQ = QQ(1:R*R,:)';
 x = VV*X;
 p = PP*P;
 non = @(x,u,p) VV*NON(UU*x,u,QQ*p);
 out = @(x,u,p) OUT(UU*x,u,QQ*p);
 OFFLINE = cputime - OFFLINE

% ONLINE
 ONLINE = cputime;
 y = rk2(non,out,T,x,U,p);
 ONLINE = cputime - ONLINE

%%%%%%%% Output %%%%%%%%

% TERMINAL
 ERROR = norm(norm(Y - y)./norm(Y))
 RELER = abs(Y - y)./abs(Y);

% PLOT
 if(nargin<1 || o==0 ), return; end
 l = (1:-0.01:0)'; cmap = [l,l,ones(101,1)];
 figure('PaperSize',[2.4,6.4],'PaperPosition',[0,0,6.4,2.4]);
 imagesc(RELER); caxis([0 max(max(RELER))]); colorbar; colormap(cmap); set(gca,'YTick',1:N);
 if(o==2 && exist('OCTAVE_VERSION')), print -dsvg nonlinear.svg; end

%%%%%%%% Integrator %%%%%%%%

function y = rk2(f,g,t,x,u,p)

 h = t(2);
 T = (t(3)-t(1))/h;
 y = zeros(numel(g(x,u(:,1),p)),T);

 for t=1:T
  x = x + h*f(x + 0.5*h*f(x,u(:,t),p),u(:,t),p); %Improved Eulers Method
  y(:,t) = g(x,u(:,t),p);
 end
