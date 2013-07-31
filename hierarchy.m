function hierarchy(o)
% hierarchical network reduction
% by Christian Himpe, 2013 ( http://gramian.de )
% released under BSD 2-Clause License ( opensource.org/licenses/BSD-2-Clause )
%*

if(exist('emgr')~=2) disp('emgr framework is required. Download at http://gramian.de/emgr.m'); return; end

%%%%%%%% Setup %%%%%%%%

 D = 3;		%Tree depth
 M = 3;		%Children per node
 S = -0.2*M;	%Diagonal Scale

 J = 1;
 O = D*M;
 N = (M^(D+1)-1)/(M-1);
 R = O;
 T = [0 0.01 1];
 L = (T(3)-T(1))/T(2);
 U = [N zeros(1,L-1)];
 X =    ones(N,1);

 A = trasm(D,M,S);
 B = sparse(N,1); B(1,1) = 1;
 C = [sparse(O,N-O) speye(O)];

 LIN = @(x,u,p) A*x + B*u;
 OUT = @(x,u,p) C*x;

%%%%%%%% Reduction %%%%%%%%

% FULL
 tic; Y = rk2(LIN,OUT,T,X,U,0); FULL = toc

% OFFLINE
 tic;
 WC = emgr(LIN,OUT,[J N O],T,'c');
 WO = emgr(LIN,OUT,[J N O],T,'o');
 [UU D VV] = squareroot(WC,WO,R);
 a = UU*A*VV;
 b = UU*B;
 c = C*VV;
 x = UU*X;
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
 if(o==2 && exist('OCTAVE_VERSION')), print -dsvg hierarchy.svg; end

%%%%%%%% Tree Assembler %%%%%%%%

function A = trasm(d,c,s)

 a = (c^(d+1)-1)/(c-1);
 A = sparse(a,a); A(1,1) = -s;

 for I=0:((a-1)/c)-1
  b = 1+c*I;
  %A(1+I,b+1:b+c) = rand(1,c)+1;
  A(b+1:b+c,1+I) = rand(c,1)+1;
  A(sub2ind(size(A),b+1,b+1):a+1:sub2ind(size(A),b+c,b+c)) = s;
 end

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

function [X Y Z] = squareroot(WC,WO,R)

 [L D l] = svd(WC); LC = L*diag(sqrt(diag(D)));
 [L D l] = svd(WO); LO = L*diag(sqrt(diag(D)));
 [U Y V] = svd(LO'*LC);
 X = ( LO*U(:,1:R)*diag(1./sqrt(diag(Y(1:R,1:R)))) )';
 Z =   LC*V(:,1:R)*diag(1./sqrt(diag(Y(1:R,1:R))));
