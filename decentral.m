function decentral(o)
% decentral (-ized control)
% by Christian Himpe, 2013-2014 (http://gramian.de )
% released under BSD 2-Clause License ( opensource.org/licenses/BSD-2-Clause )
%*

if(exist('emgr')~=2) disp('emgr framework is required. Download at http://gramian.de/emgr.m'); return; end

%%%%%%%% Setup %%%%%%%%

 J = 4;
 O = J;
 N = 32;
 T = [0 0.01 1];
 L = (T(3)-T(1))/T(2);
 U = [ones(J,1) zeros(J,L-1)];
 X =  ones(N,1);

 rand('seed',1009);
 A = rand(N,N); A(1:N+1:end) = -0.55*N; A = 0.5*(A+A');
 B = rand(N,J);
 C = B';

 LIN = @(x,u,p) A*x + B*u;
 OUT = @(x,u,p) C*x;

%%%%%%%% Reduction %%%%%%%%%

% FULL
 tic;
 Y = ab2(LIN,OUT,T,X,U,0);
 FULL = toc

% OFFLINE
 tic;
 WX = cell(J,O);
 for V=1:J
  for W=1:O
   lin = @(x,u,p) A*x + B(:,V)*u;
   out = @(x,u,p) C(W,:)*x;
   WX{V,W} = emgr(lin,out,[1 N 1],T,'x');
  end
 end

 PM = zeros(J,O);
 PM = cellfun(@(w) trace(w),WX);
 PM = PM./sum(sum(PM))

 c = zeros(J,2);
 for K=1:J
  [a b] = max(PM(:));
  [c(K,1) c(K,2)] = ind2sub(size(PM),b);
  PM(c(K,1),:) = 0;
  PM(:,c(K,2)) = 0;
 end
 OFFLINE = toc

% ONLINE
 tic;
 y = zeros(O,L);
 for K=1:J
  lin = @(x,u,p) A*x + B(:,c(K,1))*u;
  out = @(x,u,p) C(c(K,2),:)*x;
  y(K,:) = ab2(lin,out,T,X,U(c(K,1),:),0);
 end
 ONLINE = toc

%%%%%%%% Output %%%%%%%

% TERMINAL
 ERROR = norm(norm(Y - y)./norm(Y))
 RELER = abs(Y - y)./abs(Y);

% PLOT
 l = (1:-0.01:0)'; cmap = [l,l,ones(101,1)]; cmax = max(max(RELER));
 figure('PaperSize',[2.4,6.4],'PaperPosition',[0,0,6.4,2.4]);
 imagesc(RELER); caxis([0 cmax]); cbr = colorbar; colormap(cmap); 
 set(gca,'YTick',1:N,'xtick',[]); set(cbr,'YTick',[0 cmax],'YTickLabel',{'0',sprintf('%0.1e',cmax)});
 if(nargin>0), print -dsvg decentral.svg; end

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

