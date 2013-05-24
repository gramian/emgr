function advection(o)
% advection (pde transport equation reduction)
% by Christian Himpe, 2013 ( http://gramian.de )
% released under BSD 2-Clause License ( http://gramian.de/#license )
%*

if(exist('emgr')~=2) disp('emgr framework is required. Download at http://gramian.de/emgr.m'); return; end

%%%%%%%% Setup %%%%%%%%

%% CONSTANTS
 N = 100;
 J = 0;
 O = N;
 R = N/10;
 t = [0 0.01 1];
 T = (t(3)-t(1))/t(2);
 u = zeros(1,T);
 h = 0.1;
 x = 0:h:9.9; x = exp(-(x-1).^2)';
%%

%% PARAMETER
 p = 5;
 A = -eye(N) + diag(ones(N-1,1),-1); A = A*(1/h);
 C = eye(N);
%%

%% FUNCTIONS
 LIN = @(x,u,p) p*A*x;
 OUT = @(x,u,p) C*x;
%%

%%%%%%%% Linear State Reduction %%%%%%%%%

%% FULL
 %% ONLINE
  tic;
  Y = rk2(LIN,OUT,[J N O],t,x,u,p);
  ORIGINAL = toc
 %%
%%

%% WX
 %% OFFLINE
  tic;
  WO = emgr(LIN,OUT,[J N O],p,t,'o');
  [UU D VV] = svds(WO,R); VV = VV';
  a = VV*A*UU;
  c = C*UU;
  lin = @(x,u,p) p*a*x;
  out = @(x,u,p) c*x;
  OFFLINE = toc
 %%

 %% ONLINE
  tic;
  y = rk2(lin,out,[J R O],t,VV*x,u,p);
  ONLINE = toc
  ERROR = norm(norm(Y - y)./norm(Y))
 %%
%%

%% PLOT
 if(nargin<1 || o==0 ) return; end
 l = (1:-0.01:0)'; cmap = [l,l,ones(101,1)];
 figure('PaperSize',[4.8,4.8],'PaperPosition',[0,0,4.8,4.8]);
 imagesc(y); colormap(cmap);
 set(gca,'YTick',0); ylabel('X');
 if(o==2 &&  exist('OCTAVE_VERSION'))
 	print -dsvg advection.svg;
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
