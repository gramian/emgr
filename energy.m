function energy(o)
% energy (nonlinear output)
% by Christian Himpe, 2013 ( http://gramian.de )
% released under BSD 2-Clause License ( http://gramian.de/#license )
%*

if(exist('emgr')~=2) disp('emgr framework is required. Download at http://gramian.de/emgr.m'); return; end

%%%%%%%% Setup %%%%%%%%

%% CONSTANTS
 N = 10;
 t = [0 0.01 1];
 T = (t(3)-t(1))/t(2);
 u = [1 zeros(1,T-1)];
 x = N*ones(N,1);
%%

%% SYSTEM
 A = rand(N,N);
 A(1:N+1:end) = -0.48*N;
 B = rand(N,1);
 p = A(:);
%%

%% FUNCTIONS
 LIN = @(x,u,p) reshape(p,[N N])*x+B*u;
 OUT = @(x,u,p) x'*x;
%%

%%%%%%%% Parameter Reduction %%%%%%%%%

%% FULL
 %% ONLINE
  tic;
  Y = rk2(LIN,OUT,[1 N 1],t,x,u,p);
  ORIGINAL = toc
 %%
%%

%% WI
 %% OFFLINE
  tic;
  WI = emgr(LIN,OUT,[1 N 1],p,t,'i');
  [PP D QQ] = svd(WI{2}); PP = PP(1,:); QQ = QQ(1,:)'; q = PP*p;
  OFFLINE = toc
 %%

 %% ONLINE
  tic;
  y = rk2(LIN,OUT,[1 N 1],t,x,u,QQ*q);
  ONLINE = toc
  ERROR = norm(norm(Y - y)./norm(Y))
  RELER = abs(Y - y)./abs(Y);
 %%
%%

RELER = [RELER;RELER];

%% PLOT
if(nargin<1 || o==0 ) return; end
 l = (1:-0.01:0)'; cmap = [l,l,ones(101,1)];
 figure('PaperSize',[2.4,6.4],'PaperPosition',[0,0,6.4,2.4]);
 imagesc(RELER); caxis([0 max(max(RELER))]); colorbar; colormap(cmap);
 set(gca,'YTickLabel',{'','','1','',''});
 if(o==2 &&  exist('OCTAVE_VERSION'))
 	print -dsvg energy.svg;
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
