function hierarchy(o)
% hierarchy
% by Christian Himpe, 2013 ( http://gramian.de )
% released under BSD 2-Clause License ( http://gramian.de/#license )
%*

if(exist('emgr')~=2) disp('emgr framework is required. Download at http://gramian.de/emgr.m'); return; end

%%%%%%%% Setup %%%%%%%%

%% CONSTANTS
 L = 3; %Tree depth
 M = 3; %Children per node
 S = -0.2*M; %Diagonal Scale

 J = 1;
 O = L*M;
 N = (M^(L+1)-1)/(M-1);
 R = O;
 t = [0 0.01 1];
 T = (t(3)-t(1))/t(2);
 u = [N zeros(1,T-1)];
 x =    ones(N,1);
%%

%% PARAMETER
 A = trasm(L,M,S);
 B = sparse(N,1); B(1,1) = 1;
 C = [sparse(O,N-O) speye(O)];
%%

%% FUNCTIONS
 LIN = @(x,u,p) A*x + B*u;
 OUT = @(x,u,p) C*x;
%%

%%%%%%%% Hierarchy Reduction %%%%%%%%

%% FULL
 %% ONLINE
  tic;
  Y = rk2(LIN,OUT,[J N O],t,x,u,0);
  ORIGINAL = toc
 %%
%%

%% WCWO
 %% OFFLINE
 tic;
  WC = emgr(LIN,OUT,[J N O],0,t,'c');
  WO = emgr(LIN,OUT,[J N O],0,t,'o');
  [UU D VV] = squareroot(WC,WO,R);
  a = UU*A*VV;
  b = UU*B;
  c = C*VV;
  lin = @(x,u,p) a*x + b*u;
  out = @(x,u,p) c*x;
  OFFLINE = toc
 %%

 %% ONLINE
  tic;
  y = rk2(lin,out,[J R O],t,UU*x,u,0);
  ONLINE = toc
  ERROR = norm(norm(Y - y)./norm(Y))
  RELER = abs(Y - y)./abs(Y);
 %%
%%

%% PLOT
 if(nargin<1 || o==0 ) return; end
 l = (1:-0.01:0)'; cmap = [l,l,ones(101,1)];
 figure('PaperSize',[2.4,6.4],'PaperPosition',[0,0,6.4,2.4]);
 imagesc(RELER); caxis([0 max(max(RELER))]); colorbar; colormap(cmap);
 set(gca,'YTick',1:N);
 if(o==2 &&  exist('OCTAVE_VERSION'))
 	print -dsvg hierarchy.svg;
 end
%%

%%%%%%%% Tree Assembler %%%%%%%%

function A = trasm(d,c,s)

a = (c^(d+1)-1)/(c-1);
A = sparse(a,a);

A(1,1) = -s;

for I=0:((a-1)/c)-1
	b = 1+c*I;
	%A(1+I,b+1:b+c) = rand(1,c)+1;
	A(b+1:b+c,1+I) = rand(c,1)+1;
	A(sub2ind(size(A),b+1,b+1):a+1:sub2ind(size(A),b+c,b+c)) = s;
end

%%%%%%%% RK2 Integration %%%%%%%%

function y = rk2(f,g,q,t,x,u,p)
T = (t(3)-t(1))/t(2);
y = zeros(q(3),T);
h = t(2);

for A=1:T
	x = x + h*f(x + 0.5*h*f(x,u(:,A),p),u(:,A),p); %Improved Eulers Method
	y(:,A) = g(x,u(:,A),p);
end

%%%%%%%% Squareroot %%%%%%%%

function [X Y Z] = squareroot(WC,WO,R)
	[L D l] = svd(WC); LC = L*diag(sqrt(diag(D)));
	[L D l] = svd(WO); LO = L*diag(sqrt(diag(D)));
	[U Y V] = svds(LO'*LC,R);
	X = ( LO*U*diag(1./sqrt(diag(Y))) )';
	Z =   LC*V*diag(1./sqrt(diag(Y)));
