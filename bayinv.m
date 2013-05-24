function bayinv(o)
% bayinv (bayesian inverse problem combined reduction)
% by Christian Himpe, 2013 ( http://gramian.de )
% released under BSD 2-Clause License ( http://gramian.de/#license )
%*

if(exist('emgr')~=2) disp('emgr framework is required. Download at http://gramian.de/emgr.m'); return; end

warning off all
optopt = optimset('Display','off','MaxIter',10);

%%%%%%%% Setup %%%%%%%%

%% CONSTANTS
 J = 4;
 O = J;
 N = 16;
 R = O;
 t = [0 0.01 1];
 T = (t(3)-t(1))/t(2);
 u = [N*ones(J,1) zeros(J,T-1)];
 x =    ones(N,1);
%%

%% PARAMETER
 A = rand(N,N);
 A(1:N+1:end) = -0.55*N;
 A = 0.5*(A+A');
 B = rand(N,J);
 C = B';
 P = A(:);
 p = -0.5*N*eye(N);
 p = p(:);
%%

%% FUNCTIONS
 non = @(x,u,p) reshape(p,[N N])*tanh(x) + B*u;
 out = @(x,u,p) C*x;
 Y = rk2(non,out,[J N O],t,x,u,P);
%%

%%%%%%%% Nonlinear Combined Reduction %%%%%%%%

%% FULL
 %% ONLINE
  tic;
  s = @(z,w) norm2(Y-rk2(@(x,u,p) non(x,u,p),@(x,u,p) out(x,u,p),[J N O],t,x,u,z));
  r = fminunc(s,p,optopt);
  y = rk2(non,out,[J N O],t,x,u,r);
  FULL = toc
  ERROR_FULL = norm(norm(Y - y)./norm(Y))
 %%
%%

%% WJ
 %% OFFLINE
  tic;
  WJ = emgr(non,out,[J N O],p,t,'j',[0 0 0 0 0 0 0 0 0 1],1,0,x);
  [TT D VV] = svds(WJ{1},R); VV = VV';
  [PP D QQ] = svds(WJ{2},R); PP = PP';
  OFFLINE = toc
 %%

 %% ONLINE
  tic;
  s = @(z,w) norm2(Y-rk2(@(x,u,p) VV*non(TT*x,u,QQ*p),@(x,u,p) out(TT*x,u,QQ*p),[J N O],t,VV*x,u,z));
  r = fminunc(s,PP*p,optopt);
  y = rk2(@(x,u,p) VV*non(TT*x,u,p),@(x,u,p) out(TT*x,u,p),[J N O],t,VV*x,u,QQ*r);
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
	print -dsvg bayinv.svg;
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

%%%%%%%% Timeseries 2-Norm %%%%%%%%

function x = norm2(X)
	x = sqrt(sum(sum(X.*X)));
