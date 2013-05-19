function testilp(o)
% testilp (test inverse lyapunov procedure)
% by Christian Himpe, 2013 ( http://gramian.de )
% released under BSD 2-Clause License ( http://gramian.de/#license )
%*

if(exist('emgr')~=2) disp('emgr framework is required. Download at http://gramian.de/emgr.m'); return; end

%%%%%%%% Setup %%%%%%%%

%% CONSTANTS
 J = 9;
 N = 36;
 O = J;
 R = O+O;
 [A B C] = ilp(J,N,O);
 t = [0 0.01 1];
 T = (t(3)-t(1))/t(2);
 u = [ones(J,1) zeros(J,T-1)];
 x = zeros(N,1);
%%

%% FUNCTIONS
 LIN = @(x,u,p) A*x + B*u;
 OUT = @(x,u,p) C*x;
%%

%%%%%%%% Reduction %%%%%%%%%

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
  [UU D VV] = balance(WC,WO); UU = UU(1:R,:); VV = VV(:,1:R);
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
 	print -dsvg ilp.svg;
 end
%%

%%%%%%%% RK2 Integration %%%%%%%%

function y = rk2(f,g,q,t,x,u,p)
T = (t(3)-t(1))/t(2);
y = zeros(q(3),T);
h = t(2);

for A=1:T
	x = x + h*f(x + 0.5*h*f(x,u(:,A),p),u(:,A),p); %Improved Eulers Method
	y(:,A) = g(x,u(:,A),p);
end

%%%%%%%% Balance %%%%%%%%

function [X Y Z] = balance(WC,WO)
	L = chol(WC+eye(size(WC,1)))-eye(size(WC,1));
	[U Y V] = svd(L*WO*L');
	X = diag(sqrt(diag(Y))) * V' / L';
	Z = L'*U*diag(1./sqrt(diag(Y)));


%%%%%%%% Inverse Lyapunov Procedure %%%%%%%%%

function [A B C] = ilp(J,N,O,s)

%% Gramian Eigenvalues
 WC = exp(-N + N*rand(N,1));
 WO = exp(-N + N*rand(N,1));

%% Gramian Eigenvectors
 X = randn(N,N);
 [U E V] = svd(X);

%% Balancing Trafo
 [P D Q] = svd(diag(WC.*WO));
 W = -D;

%% Input and Output
 B = randn(N,J);

 if(nargin<4 || s==0)
        C = randn(O,N);
 else
        C = B';
 end

%% Scale Output Matrix
 BB = sum(B.*B,2);  % = diag(B*B')
 CC = sum(C.*C,1)'; % = diag(C'*C)
 C = bsxfun(@times,C,sqrt(BB./CC)');

%% Solve System Matrix
 f = @(x,u,p) W*x+B*u;
 g = @(x,u,p) C*x;
 A = -emgr(f,g,[J N O],0,[0 0.01 1],'c');

%% Unbalance System
 T = U'*P';
 A = T*A*T';
 B = T*B;
 C = C*T';
