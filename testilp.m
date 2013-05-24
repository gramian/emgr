function testilp(o)
% testilp (test inverse lyapunov procedure)
% by Christian Himpe, 2013 ( http://gramian.de )
% released under BSD 2-Clause License ( http://gramian.de/#license )
%*

if(exist('emgr')~=2) disp('emgr framework is required. Download at http://gramian.de/emgr.m'); return; end
if(exist('ilp')~=2)  disp('ilp generator is required. Download at http://gramian.de/ilp.m'); return; end

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
  [UU D VV] = balance(WC,WO,R);
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

function [X Y Z] = balance(WC,WO,R)
	L = chol(WC+eye(size(WC,1)))-eye(size(WC,1));
	[U Y V] = svds(L*WO*L',R);
	X = diag(sqrt(diag(Y))) * V' / L';
	Z = L'*U*diag(1./sqrt(diag(Y)));
