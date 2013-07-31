function largescale(o)
% largescale (linear state reduction)
% by Christian Himpe, 2013 ( http://gramian.de )
% released under BSD 2-Clause License ( opensource.org/licenses/BSD-2-Clause )
%*

if(exist('emgr')~=2) disp('emgr framework is required. Download at http://gramian.de/emgr.m'); return; end

%%%%%%%% Setup %%%%%%%%

 J = 16;
 O = J;
 N = 128;
 R = O;
 T = [0 0.01 1];
 L = (T(3)-T(1))/T(2);
 U = [N*ones(J,1) zeros(J,L-1)];
 X =    ones(N,1);

 rand('seed',1009);
 A = rand(N,N); A(1:N+1:end) = -0.55*N; A = 0.5*(A+A');
 B = rand(N,J);
 C = B';

 LIN = @(x,u,p) A*x + B*u;
 OUT = @(x,u,p) C*x;

%%%%%%%% Reduction %%%%%%%%%

% FULL
 tic; Y = rk2(LIN,OUT,T,X,U,0); FULL = toc

% OFFLINE
 tic;
 WX = emgr(LIN,OUT,[J N O],T,'x');
 [UU D VV] = svd(WX); UU = UU(:,1:R); VV = VV(:,1:R)';
 a = VV*A*UU;
 b = VV*B;
 c = C*UU;
 x = VV*X;
 lin = @(x,u,p) a*x + b*u;
 out = @(x,u,p) c*x;
 OFFLINE = toc

% ONLINE
 tic; y = rk2(lin,out,T,x,U,0); ONLINE = toc

%%%%%%%% Output %%%%%%%%

% TERMINAL
 ERROR = norm(norm(Y - y)./norm(Y))
 RELER = abs(Y - y)./abs(Y);

% PLOT 
 if(nargin<1 || o==0 ), return; end
 l = (1:-0.01:0)'; cmap = [l,l,ones(101,1)];
 figure('PaperSize',[2.4,6.4],'PaperPosition',[0,0,6.4,2.4]);
 imagesc(RELER); caxis([0 max(max(RELER))]); colorbar; colormap(cmap); set(gca,'YTick',1:N);
 if(o==2 && exist('OCTAVE_VERSION')), print -dsvg largescale.svg; end

%%%%%%%% Integrator %%%%%%%%

function y = rk2(f,g,t,x,u,p)

 h = t(2);
 T = (t(3)-t(1))/h;
 y = zeros(numel(g(x,u(:,1),p)),T);

 for t=1:T
  x = x + h*f(x + 0.5*h*f(x,u(:,t),p),u(:,t),p); %Improved Eulers Method
  y(:,t) = g(x,u(:,t),p);
 end
