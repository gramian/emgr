function closed(o)
% closed loop reduction
% by Christian Himpe, 2013 ( http://gramian.de )
% released under BSD 2-Clause License ( http://gramian.de/#license )
%*

if(exist('emgr')~=2) disp('emgr framework is required. Download at http://gramian.de/emgr.m'); return; end

%%%%%%%% Setup %%%%%%%%

%% CONSTANTS
 J = 8;
 O = J;
 N = 32;
 R = O;
 t = [0 0.01 1];
 T = (t(3)-t(1))/t(2);
 u = [N*ones(J,1) zeros(J,T-1)];
 x = zeros(N,1);
%%

%% PARAMETER
 A = rand(N,N);
 A(1:N+1:end) = -0.55*N;
 A = 0.5*(A+A');
 B = rand(N,J);
 C = B';
%%

%% FUNCTIONS
 NON = @(x,u,p) A*tanh(x) + B*(u+C*x);
 OUT = @(x,u,p) C*x;
%%

%%%%%%%% Nonlinear State Reduction %%%%%%%%%

%% FULL
 %% ONLINE
  tic;
  Y = rk2(NON,OUT,[J N O],t,x,u,0);
  ORIGINAL = toc
 %%
%%

%% WX
 %% OFFLINE
  tic;
  WX = emgr(NON,OUT,[J N O],0,t,'x');
  [UU D VV] = svds(WX,R); VV = VV';
  a = VV*A*UU;
  b = VV*B;
  c = C*UU;
  non = @(x,u,p) a*tanh(x) + b*(u+c*x);
  out = @(x,u,p) c*x;
  OFFLINE = toc
 %%

 %% ONLINE
  tic;
  y = rk2(non,out,[J R O],t,VV*x,u,0);
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
 	print -dsvg closed.svg;
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
 
