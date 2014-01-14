function closed(o)
% closed loop reduction
% by Christian Himpe, 2013,2014 ( http://gramian.de )
% released under BSD 2-Clause License ( opensource.org/licenses/BSD-2-Clause )
%*

if(exist('emgr')~=2) disp('emgr framework is required. Download at http://gramian.de/emgr.m'); return; end

%%%%%%%% Setup %%%%%%%%

 J = 8;
 O = J;
 N = 32;
 R = O;
 T = [0 0.01 1];
 L = (T(3)-T(1))/T(2);
 U = [N*ones(J,1) zeros(J,L-1)];
 X =   zeros(N,1);

 rand('seed',1009);
 A = rand(N,N); A(1:N+1:end) = -0.55*N; A = 0.5*(A+A');
 B = rand(N,J);
 C = B';

 NON = @(x,u,p) A*tanh(x) + B*(u+C*x);
 OUT = @(x,u,p) C*x;

%%%%%%%% Reduction %%%%%%%%%

% FULL
 FULL = cputime;
 Y = rk2(NON,OUT,T,X,U,0);
 FULL = cputime - FULL

% OFFLINE
 OFFLINE = cputime;
 WX = emgr(NON,OUT,[J N O],T,'x');
 [UU D VV] = svd(WX); UU = UU(:,1:R); VV = VV(:,1:R)';
 a = VV*A*UU;
 b = VV*B;
 c = C*UU;
 x = VV*X;
 non = @(x,u,p) a*tanh(x) + b*(u+c*x);
 out = @(x,u,p) c*x;
 OFFLINE = cputime - OFFLINE

% ONLINE
 ONLINE = cputime;
 y = rk2(non,out,T,x,U,0);
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
 if(o==2 && exist('OCTAVE_VERSION')), print -dsvg closed.svg; end

%%%%%%%% Integrator %%%%%%%%

function y = rk2(f,g,t,x,u,p)

 h = t(2);
 T = (t(3)-t(1))/h;
 y = zeros(numel(g(x,u(:,1),p)),T);

 for t=1:T
  x = x + h*f(x + 0.5*h*f(x,u(:,t),p),u(:,t),p); %Improved Eulers Method
  y(:,t) = g(x,u(:,t),p);
 end
