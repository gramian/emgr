function testilp(o)
% testilp (test inverse lyapunov procedure)
% by Christian Himpe, 2013,2014 ( http://gramian.de )
% released under BSD 2-Clause License ( opensource.org/licenses/BSD-2-Clause )
%*

if(exist('emgr')~=2) disp('emgr framework is required. Download at http://gramian.de/emgr.m'); return; end
if(exist('ilp')~=2)  disp('ilp generator is required. Download at http://gramian.de/ilp.m'); return; end

%%%%%%%% Setup %%%%%%%%

 J = 8;
 N = 64;
 O = J;
 R = O+O;
 [A B C] = ilp(J,N,O,0,1009);
 T = [0 0.01 1];
 L = (T(3)-T(1))/T(2);
 U = [ones(J,1) zeros(J,L-1)];
 X = zeros(N,1);

 LIN = @(x,u,p) A*x + B*u;
 OUT = @(x,u,p) C*x;

%%%%%%%% Reduction %%%%%%%%%

% FULL
 FULL = cputime;
 Y = rk2(LIN,OUT,T,X,U,0);
 FULL = cputime - FULL

% OFFLINE
 OFFLINE = cputime;
 WC = emgr(LIN,OUT,[J N O],T,'c');
 WO = emgr(LIN,OUT,[J N O],T,'o');
 [UU D VV] = balance(WC,WO,R);
 a = UU*A*VV;
 b = UU*B;
 c = C*VV;
 x = UU*X;
 lin = @(x,u,p) a*x + b*u;
 out = @(x,u,p) c*x;
 OFFLINE = cputime - OFFLINE

% ONLINE
 ONLINE = cputime;
 y = rk2(lin,out,T,x,U,0);
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
 if(o==2 && exist('OCTAVE_VERSION')), print -dsvg testilp.svg; end

%%%%%%%% Integrator %%%%%%%%

function y = rk2(f,g,t,x,u,p)

 h = t(2);
 T = (t(3)-t(1))/h;
 y = zeros(numel(g(x,u(:,1),p)),T);

 for t=1:T
  x = x + h*f(x + 0.5*h*f(x,u(:,t),p),u(:,t),p); %Improved Eulers Method
  y(:,t) = g(x,u(:,t),p);
 end

%%%%%%%% Balancer %%%%%%%%

function [X Y Z] = balance(WC,WO,R)

 L = chol(WC+eye(size(WC,1)))-eye(size(WC,1));
 [U Y V] = svd(L*WO*L');
 X = diag(sqrt(diag(Y(1:R,1:R)))) * V(:,1:R)' / L';
 Z = L'*U(:,1:R)*diag(1./sqrt(diag(Y(1:R,1:R))));
