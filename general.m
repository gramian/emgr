function general(o)
% general (generalized cross gramian via controllability gramian)
% by Christian Himpe, 2013 ( http://gramian.de )
% released under BSD 2-Clause License ( http://gramian.de/#license )
%*

if(exist('emgr')~=2) disp('emgr framework is required. Download at http://gramian.de/emgr.m'); return; end

%%%%%%%% Setup %%%%%%%%

%% CONSTANTS
 J = 16;
 N = 128;
 O = J;
 R = O;
 t = [0 0.01 1];
 T = (t(3)-t(1))/t(2);
 u = [N*ones(J,1) zeros(J,T-1)];
 x =  ones(N,1);
%%

%% PARAMETER
 A = rand(N,N);
 A(1:N+1:end) = -0.55*N;
 A = 0.5*(A+A');
 B = rand(N,J);
 C = B';
%%

%% FUNCTIONS
 LIN = @(x,u,p) A*x+B*u;
 OUT = @(x,u,p) C*x;

 Lin = @(x,u,p) [A,sparse(N,N);sparse(N,N),A']*x + [B;C']*u;
 Out = @(x,u,p) x;
%%

%%%%%%%% Linear State Reduction %%%%%%%%%

%% FULL
 %% ONLINE
  tic;
  Y = rk2(LIN,OUT,[J N O],t,x,u,0);
  ORIGINAL = toc
 %%
%%

%% WX
 %% OFFLINE
  tic;
  WG = emgr(Lin,Out,[J 2*N 2*N],0,t,'c');
  %WC = WG(1:N,1:N);
  %WO = WG(N+1:N+N,N+1:N+N);
  WX = WG(1:N,N+1:N+N);

  [UU D VV] = svd(WX); UU = UU(:,1:R); VV = VV(:,1:R)';
  a = VV*A*UU;
  b = VV*B;
  c = C*UU;
  lin = @(x,u,p) a*x + b*u;
  out = @(x,u,p) c*x;
  OFFLINE = toc
 %%

 %% ONLINE
  tic;
  y = rk2(lin,out,[J R O],t,VV*x,u,0);
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
 	print -dsvg general.svg;
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
 
