function energy(o)
% energy reduction (nonlinear output)
% by Christian Himpe, 2013,2014 ( http://gramian.de )
% released under BSD 2-Clause License ( opensource.org/licenses/BSD-2-Clause )
%*

if(exist('emgr')~=2) disp('emgr framework is required. Download at http://gramian.de/emgr.m'); return; end

%%%%%%%% Setup %%%%%%%%

 N = 8;
 T = [0 0.01 1];
 L = (T(3)-T(1))/T(2);
 U = [1 zeros(1,L-1)];
 X = N*ones(N,1);

 rand('seed',1009);
 A = rand(N,N); A(1:N+1:end) = -0.48*N;
 B = rand(N,1);
 P = A(:);

 LIN = @(x,u,p) reshape(p,[N N])*x+B*u;
 OUT = @(x,u,p) x'*x;

%%%%%%%% Reduction %%%%%%%%%

% FULL
 FULL = cputime;
 Y = rk2(LIN,OUT,T,X,U,P);
 FULL = cputime - FULL

% OFFLINE
 OFFLINE = cputime;
 WI = emgr(LIN,OUT,[1 N 1],T,'i',P);
 [PP D QQ] = svd(WI{2}); PP = PP(1,:); QQ = QQ(1,:)';
 p = QQ*PP*P;
 OFFLINE = cputime - OFFLINE

% ONLINE
 ONLINE = cputime;
 y = rk2(LIN,OUT,T,X,U,p);
 ONLINE = cputime - ONLINE

%%%%%%%% Output %%%%%%%%

% TERMINAL
 ERROR = norm(norm(Y - y)./norm(Y))
 RELER = abs(Y - y)./abs(Y); RELER = [RELER;RELER];

% PLOT
 if(nargin<1 || o==0 ), return; end
 l = (1:-0.01:0)'; cmap = [l,l,ones(101,1)];
 figure('PaperSize',[2.4,6.4],'PaperPosition',[0,0,6.4,2.4]);
 imagesc(RELER); caxis([0 max(max(RELER))]); colorbar; colormap(cmap); set(gca,'YTickLabel',{'','','1','',''});
 if(o==2 && exist('OCTAVE_VERSION')), print -dsvg energy.svg; end

%%%%%%%% Integrator %%%%%%%%

function y = rk2(f,g,t,x,u,p)

 h = t(2);
 T = (t(3)-t(1))/h;
 y = zeros(numel(g(x,u(:,1),p)),T);

 for t=1:T
  x = x + h*f(x + 0.5*h*f(x,u(:,t),p),u(:,t),p); %Improved Eulers Method
  y(:,t) = g(x,u(:,t),p);
 end
