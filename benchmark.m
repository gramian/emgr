function benchmark(o)
% bench (iss model reduction benchmark)
% by Christian Himpe, 2013 ( http://gramian.de )
% released under BSD 2-Clause License ( http://gramian.de/#license )
%*

if(exist('emgr')~=2) disp('emgr framework is required. Download at http://gramian.de/emgr.m'); return; end

%%%%%%%% Setup %%%%%%%%

%% CONSTANTS
 D = 'iss';
 if(exist(['/tmp/',D,'.mat'],'file')==0) unzip(['http://www.slicot.org/shared/bench-data/',D,'.zip'],'/tmp'); end
 load(['/tmp/',D,'.mat']);
 N = 270;
 n = N/2;
 J = 3;
 O = 3;
 R = 26;
 r = R/2;
 t = [0 0.001 1];
 T = (t(3)-t(1))/t(2);
 u = [ones(J,1) zeros(J,T-1)];
 x = zeros(N,1);
%%

%% FUNCTIONS
 LIN = @(x,u,p) A*x + B*u;
 OUT = @(x,u,p) C*x;
%%

%%%%%%%% Linear State Reduction %%%%%%%%%

%% FULL
 %% ONLINE
  tic;
  Y = rk2(LIN,OUT,[J N O],[0 0.01 1],x,u,0);
  ORIGINAL = toc
 %%
%%

%% WX
 %% OFFLINE
  tic;
  WX = emgr(LIN,OUT,[J N O],0,t,'x',[0 0 0 0 0 0 0 0 0 1]);
  %WXP = WX(1:n,1:n);
  %[U D V] = svds(WXP,r); V = U';
  WXV = WX(n+1:N,n+1:N);
  [U D V] = svds(WXV,r); V = U';

  a = [zeros(r,r),eye(r);V*A(n+1:N,1:n)*U,V*A(n+1:N,n+1:N)*U];
  b = [zeros(r,J);V*B(n+1:N,:)];
  c = [zeros(O,r),C(:,n+1:N)*U];

  lin = @(x,u,p) a*x + b*u;
  out = @(x,u,p) c*x;
  OFFLINE = toc
 %%

 %% ONLINE
  tic;
  y = rk2(lin,out,[J R O],[0 0.01 1],zeros(R,1),u,0);
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
 	print -dsvg benchmark.svg;
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
