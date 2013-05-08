function goalor(o)
% goalor
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
 B = rand(N,J);
 C = rand(O,N);
 P = A(:);
 p = -0.5*N*eye(N);
 p = p(:);
%%

%% FUNCTIONS
 LIN = @(x,u,p) reshape(p,[N N])*x + B*u;
 OUT = @(x,u,p) C*x;
 Y = rk2(LIN,OUT,[J N O],t,x,u,P);
 YD = cell(2,1); YD{2,1} = Y;
 Z = cell(2,1);
%%

%%%%%%%% Goal-Oriented Inference %%%%%%%%

%% FULL
 %% ONLINE
  tic;
  s = @(z,w) norm2(Y-rk2(LIN,OUT,[J N O],t,x,u,z));
  r = fminunc(s,p,optopt);
  y = rk2(LIN,OUT,[J N O],t,x,u,r);
  FULL_offline_online_error = [0,toc,norm(norm(Y - y)./norm(Y))]
  Z{1} = abs(Y - y)./abs(Y);
 %%
%%

%% BOWO
 %% OFFLINE
  tic;
  BO = emgr(LIN,OUT,[J N O],p,t,'o',[0 0 0 3 0 3 0 0 1 0],1,0,x,1,1,YD);
  WO = emgr(LIN,OUT,[J N O],p,t,'o',0,0,0,x);
  [UU D VV] = balance(BO,WO); UU = UU(1:R,:); VV = VV(:,1:R);
  b = UU*B;
  c = C*VV;
  q = -0.5*R*eye(R);
  q = q(:);
  lin = @(x,u,p) reshape(p,[R R])*x + b*u;
  out = @(x,u,p) c*x;
  OFFLINE = toc;
 %%

 %% ONLINE
  tic;
  s = @(z,w) norm2(Y-rk2(lin,out,[J R O],t,UU*x,u,z));
  r = fminunc(s,q,optopt);
  y = rk2(lin,out,[J R O],t,UU*x,u,r);
  BOWO_offline_online_error = [OFFLINE,toc,norm(norm(Y - y)./norm(Y))]
  Z{2} = abs(Y - y)./abs(Y);
 %%
%%

%% PLOT
 if(nargin<1 || o==0 ) return; end
 l = (1:-0.01:0)'; cmap = [l,l,ones(101,1)];
 m = max([max(max(Z{1})),max(max(Z{2}))]);

 f0 = figure('PaperSize',[1.6,6.4],'PaperPosition',[0,0,6.4,1.6]);
 imagesc(Z{1}); xlabel('Time'); ylabel('Relative Error'); caxis([0 m]); colorbar; colormap(cmap);
 if(o==2) saveas(f0,'base','epsc'); end

 f3 = figure('PaperSize',[1.6,6.4],'PaperPosition',[0,0,6.4,1.6]);
 imagesc(Z{2}); xlabel('Time'); ylabel('Relative Error'); caxis([0 m]); colorbar; colormap(cmap);
 if(o==2) saveas(f3,'bowo','epsc'); end
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
	[U D V] = svd(L*WO*L');
	X = diag(diag(D).^0.25) * V' * inv(L');
	Y = X'*WO*X;
	Z = inv(X);

%%%%%%%% Timeseries 2-Norm %%%%%%%%

function x = norm2(X)
	x = sqrt(sum(sum(X.*X)));
