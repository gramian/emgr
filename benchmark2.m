function benchmark2(o)
% nonlinear benchmark (nonlinear rc ladder)
% by Christian Himpe, 2013 ( http://gramian.de )
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
 U = [zeros(1,4) ones(1,L-4)];
 X = zeros(N,1);

 g = @(x) exp(x)+x-1.0;
 
 A1 = spdiags(ones(N-1,1),-1,N,N)-speye(N);
 A2 = spdiags([ones(N-1,1);0],0,N,N)-spdiags(ones(N,1),1,N,N);
 
 NON = @(x,u,p) g(A1*x)-g(A2*x) + [u;sparse(N-1,1)];
 OUT = @(x,u,p) x(1);

%%%%%%%% Reduction %%%%%%%%%

% FULL
 tic; Y = rk2(NON,OUT,T,X,U,0); FULL = toc

% OFFLINE
 tic;
 WX = emgr(NON,OUT,[J N O],T,'x')
 [UU D VV] = svd(WX); UU = UU(:,1:R); VV = VV(:,1:R)'; diag(D)
 x = VV*X;
 non = @(x,u,p) VV*NON(UU*x,u,p);
 out = @(x,u,p) OUT(UU*x,u,p);
 OFFLINE = toc

% ONLINE
 tic; y = rk2(non,out,T,x,U,0); ONLINE = toc

%%%%%%%% Output %%%%%%%%

% TERMINAL
 ERROR = norm(norm(Y - y)./norm(Y))
 RELER = abs(Y - y)./abs(Y); RELER = [RELER;RELER];

% PLOT
 if(nargin<1 || o==0 ), return; end
 l = (1:-0.01:0)'; cmap = [l,l,ones(101,1)];
 figure('PaperSize',[2.4,6.4],'PaperPosition',[0,0,6.4,2.4]);
 imagesc(RELER); caxis([0 max(max(RELER))]); colorbar; colormap(cmap); set(gca,'YTickLabel',{'','','1','',''});
 if(o==2 && exist('OCTAVE_VERSION')), print -dsvg benchmark2.svg; end

%%%%%%%% Integrator %%%%%%%%

function y = rk2(f,g,t,x,u,p)

 h = t(2);
 T = (t(3)-t(1))/h;
 y = zeros(numel(g(x,u(:,1),p)),T);

 for t=1:T
  x = x + h*f(x + 0.5*h*f(x,u(:,t),p),u(:,t),p); %Improved Eulers Method
  y(:,t) = g(x,u(:,t),p);
 end
