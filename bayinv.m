function bayinv(o)
% bayinv (bayesian inverse problem combined reduction)
% by Christian Himpe, 2013 ( http://gramian.de )
% released under BSD 2-Clause License ( opensource.org/licenses/BSD-2-Clause )
%*

if(exist('emgr')~=2) disp('emgr framework is required. Download at http://gramian.de/emgr.m'); return; end

%%%%%%%% Setup %%%%%%%%

 J = 4;
 O = J;
 N = 16;
 R = O;
 T = [0 0.01 1];
 L = (T(3)-T(1))/T(2);
 U = [N*ones(J,1) zeros(J,L-1)];
 X =    ones(N,1);

 rand('seed',1009);
 A = rand(N,N); A(1:N+1:end) = -0.55*N; A = 0.5*(A+A');
 B = rand(N,J);
 C = B';
 P = A(:);
 Q = -0.5*N*eye(N); Q = Q(:);

 NON = @(x,u,p) reshape(p,[N N])*tanh(x) + B*u;
 OUT = @(x,u,p) C*x;

 warning off all
 OPT = optimset('Display','off','MaxIter',10);

%%%%%%%% Reduction %%%%%%%%

% ORIGINAL
 Y = rk2(NON,OUT,[J N O],T,X,U,P);

% FULL
 tic;
 j = @(q,o) norm(Y-rk2(NON,OUT,[J N O],T,X,U,q));
 r = fminunc(j,Q,OPT);
 y = rk2(NON,OUT,[J N O],T,X,U,r); FULL = toc
 ERROR_FULL = norm(norm(Y - y)./norm(Y))

% OFFLINE
 tic;
 WJ = emgr(NON,OUT,[J N O],Q,T,'j',[0 0 0 0 0 0 0 0 0 1],1,0,X);
 [TT D VV] = svd(WJ{1}); TT = TT(:,1:R);   VV = VV(:,1:R)';
 [PP D QQ] = svd(WJ{2}); PP = PP(1:R*R,:); QQ = QQ(1:R*R,:)';
 x = VV*X;
 q = PP*Q;
 non = @(x,u,p) VV*NON(TT*x,u,QQ*p);
 out = @(x,u,p) OUT(TT*x,u,QQ*p);
 OFFLINE = toc

% ONLINE
 tic;
 j = @(q,o) norm(Y-rk2(non,out,[J N O],T,x,U,q));
 r = fminunc(j,q,OPT);
 y = rk2(non,out,[J N O],T,x,U,r); ONLINE = toc

%%%%%%%% Output %%%%%%%%

% TERMINAL
 ERROR = norm(norm(Y - y)./norm(Y))
 RELER = abs(Y - y)./abs(Y);

% PLOT
 if(nargin<1 || o==0 ), return; end
 l = (1:-0.01:0)'; cmap = [l,l,ones(101,1)];
 figure('PaperSize',[2.4,6.4],'PaperPosition',[0,0,6.4,2.4]);
 imagesc(RELER); caxis([0 max(max(RELER))]); colorbar; colormap(cmap); set(gca,'YTick',1:N);
 if(o==2 && exist('OCTAVE_VERSION')), print -dsvg bayinv.svg; end

%%%%%%%% Integrator %%%%%%%%

function y = rk2(f,g,q,t,x,u,p)

 T = (t(3)-t(1))/t(2);
 y = zeros(q(3),T);
 h = t(2);

 for A=1:T
  x = x + h*f(x + 0.5*h*f(x,u(:,A),p),u(:,A),p); %Improved Eulers Method
  y(:,A) = g(x,u(:,A),p);
 end
