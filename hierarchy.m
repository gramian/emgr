function hierarchy(o)
% hierarchy (hierarchical network reduction)
% by Christian Himpe, 2013-2014 ( http://gramian.de )
% released under BSD 2-Clause License ( opensource.org/licenses/BSD-2-Clause )
%*

if(exist('emgr')~=2) disp('emgr framework is required. Download at http://gramian.de/emgr.m'); return; end

%%%%%%%% Setup %%%%%%%%

 D = 3;		%Tree depth
 M = 2;		%Children per node

 J = 1;
 O = M^D;
 N = (M^(D+1)-1)/(M-1);
 R = M^(D-1);
 T = [0 1 100];
 L = (T(3)-T(1))/T(2);
 U = exp(-0.0005*(T(2):T(3)).^2);
 X = zeros(N,1);

 rand('seed',1009);
 A = trasm(D,M);
 B = sparse(N,1); B(1,1) = D;
 C = [sparse(O,N-O) speye(O)];

 LIN = @(x,u,p) A*x + B*u;
 OUT = @(x,u,p) C*x;

 ADJ = @(x,u,p) A'*x + C'*u;
 AOU = @(x,u,p) B'*x;

%%%%%%%% Reduction %%%%%%%%

% FULL
 tic;
 Y = rk1(LIN,OUT,T,X,U,0);
 FULL = toc

% OFFLINE
 tic;
 WC = emgr(LIN,OUT,[J N O],T,'c');
 WO = emgr(ADJ,AOU,[O N J],T,'c');
 [UU D VV] = squareroot(WC,WO,R);
 a = UU*A*VV;
 b = UU*B;
 c = C*VV;
 x = UU*X;
 lin = @(x,u,p) a*x + b*u;
 out = @(x,u,p) c*x;
 OFFLINE = toc

% ONLINE
 tic;
 y = rk1(lin,out,T,x,U,0);
 ONLINE = toc

%%%%%%%% Output %%%%%%%%

% TERMINAL
 ERROR = norm(norm(Y - y)./norm(Y))
 RELER = abs(Y - y);

% PLOT
 if(nargin==0), return; end
 l = (1:-0.01:0)'; cmap = [l,l,ones(101,1)]; cmax = max(max(RELER));
 figure('PaperSize',[2.4,6.4],'PaperPosition',[0,0,6.4,2.4]);
 imagesc(RELER); caxis([0 cmax]); cbr = colorbar; colormap(cmap); 
 set(gca,'YTick',1:N,'xtick',[]); set(cbr,'YTick',[0 cmax],'YTickLabel',{'0',sprintf('%0.1e',cmax)});
 print -dsvg hierarchy.svg;

%%%%%%%% Tree Assembler %%%%%%%%

function A = trasm(d,c)

 a = (c^(d+1)-1)/(c-1);
 A = -speye(a);

 for I=0:((a-1)/c)-1
  b = 1+c*I;
  %A(1+I,b+1:b+c) = rand(1,c)+1;
  A(b+1:b+c,1+I) = rand(c,1)+1;
 end

%%%%%%%% Integrator %%%%%%%%

function y = rk1(f,g,t,x,u,p)

 h = t(2);
 T = (t(3)-t(1))/h;
 y = zeros(numel(g(x,u(:,1),p)),T);

 for t=1:T
  x = x + h*f(x,u(:,t),p);
  y(:,t) = g(x,u(:,t),p);
 end

%%%%%%%% Balancer %%%%%%%%

function [X Y Z] = squareroot(WC,WO,R)

 [L D l] = svd(WC); LC = L*diag(sqrt(diag(D)));
 [L D l] = svd(WO); LO = L*diag(sqrt(diag(D)));
 [U Y V] = svd(LO'*LC);
 X = ( LO*U(:,1:R)*diag(1.0./sqrt(diag(Y(1:R,1:R)))) )';
 Z =   LC*V(:,1:R)*diag(1.0./sqrt(diag(Y(1:R,1:R))));
