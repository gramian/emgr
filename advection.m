function advection(o)
% pde transport equation reduction
% by Christian Himpe, 2013,2014 ( http://gramian.de )
% released under BSD 2-Clause License ( opensource.org/licenses/BSD-2-Clause )
%*

if(exist('emgr')~=2) disp('emgr framework is required. Download at http://gramian.de/emgr.m'); return; end

%%%%%%%% Setup %%%%%%%%

 N = 100;
 J = 0;
 O = N;
 R = N/10;
 T = [0 0.01 1];
 L = (T(3)-T(1))/T(2);
 U = zeros(1,L);
 H = 0.1;
 X = 0:H:9.9; X = exp(-(X-1).^2)';

 P = 5;
 A = spdiags((1.0/H)*[ones(N,1) -ones(N,1)],[-1,0],N,N);

 LIN = @(x,u,p) p*A*x;
 OUT = @(x,u,p) x;

%%%%%%%% Reduction %%%%%%%%%

% FULL
 FULL = cputime;
 Y = rk2(LIN,OUT,T,X,U,P);
 FULL = cputime - FULL

% OFFLINE
 OFFLINE = cputime;
 WO = emgr(LIN,1,[J N O],T,'o',P);
 [UU D VV] = svd(WO); UU = UU(:,1:R); VV = VV(:,1:R)';
 a = VV*A*UU;
 c = UU;
 x = VV*X;
 lin = @(x,u,p) p*a*x;
 out = @(x,u,p) c*x;
 OFFLINE = cputime - OFFLINE

% ONLINE
 ONLINE = cputime;
 y = rk2(lin,out,T,x,U,P);
 ONLINE = cputime - ONLINE

%%%%%%%% Output %%%%%%%

% TERMINAL
 ERROR = norm(norm(Y - y)./norm(Y))

% PLOT
 if(nargin<1 || o==0 ), return; end
 l = (1:-0.01:0)'; cmap = [l,l,ones(101,1)];
 figure('PaperSize',[4.8,4.8],'PaperPosition',[0,0,4.8,4.8]);
 imagesc(y); colormap(cmap); set(gca,'YTick',0); ylabel('X');
 if(o==2 && exist('OCTAVE_VERSION')), print -dsvg advection.svg; end

%%%%%%%% Integrator %%%%%%%%

function y = rk2(f,g,t,x,u,p)

 h = t(2);
 T = (t(3)-t(1))/h;
 y = zeros(numel(g(x,u(:,1),p)),T);

 for t=1:T
  x = x + h*f(x + 0.5*h*f(x,u(:,t),p),u(:,t),p); %Improved Eulers Method
  y(:,t) = g(x,u(:,t),p);
 end
