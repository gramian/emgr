function benchmark(o)
% bench (iss model reduction benchmark)
% by Christian Himpe, 2013 ( http://gramian.de )
% released under BSD 2-Clause License ( opensource.org/licenses/BSD-2-Clause )
%*

if(exist('emgr')~=2) disp('emgr framework is required. Download at http://gramian.de/emgr.m'); return; end

%%%%%%%% Download %%%%%%%%

 D = 'iss';
 if(exist(['/tmp/',D,'.mat'],'file')==0) unzip(['http://www.slicot.org/shared/bench-data/',D,'.zip'],'/tmp'); end
 load(['/tmp/',D,'.mat']);

%%%%%%%% Setup %%%%%%%%

 N = 270;
 n = N/2;
 J = 3;
 O = 3;
 R = 26;
 r = R/2;
 T = [0 0.01 1];
 L = (T(3)-T(1))/T(2);
 U = [ones(J,1) zeros(J,L-1)];
 X = zeros(N,1);

 LIN = @(x,u,p) A*x + B*u;
 OUT = @(x,u,p) C*x;

%%%%%%%% Reduction %%%%%%%%%

% FULL
 tic; Y = rk2(LIN,OUT,T,X,U,0); FULL = toc

% OFFLINE
 tic;
 WX = emgr(LIN,OUT,[J N O],T,'x',0,[0 0 0 0 0 0 0 0 0 2]);
 %WXP = WX(1:n,1:n);
 %[UU D VV] = svd(WXP); UU = UU(:,1:r); VV = UU';
 WXV = WX(n+1:N,n+1:N);
 [UU D VV] = svd(WXV); UU = UU(:,1:r); VV = UU';
 a = [zeros(r,r),eye(r);VV*A(n+1:N,1:n)*UU,VV*A(n+1:N,n+1:N)*UU];
 b = [zeros(r,J);VV*B(n+1:N,:)];
 c = [zeros(O,r),C(:,n+1:N)*UU];
 x = zeros(R,1);
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
 if(o==2 && exist('OCTAVE_VERSION')), print -dsvg benchmark.svg; end

%%%%%%%% Integrator %%%%%%%%

function y = rk2(f,g,t,x,u,p)

 h = t(2);
 T = (t(3)-t(1))/h;
 y = zeros(numel(g(x,u(:,1),p)),T);

 for t=1:T
  x = x + h*f(x + 0.5*h*f(x,u(:,t),p),u(:,t),p); %Improved Eulers Method
  y(:,t) = g(x,u(:,t),p);
 end
