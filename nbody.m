function nbody(s)
% nbody (n-body reduction)
% by Christian Himpe, 2013-2014 ( http://gramian.de )
% released under BSD 2-Clause License ( opensource.org/licenses/BSD-2-Clause )
%*

if(exist('emgr')~=2) disp('emgr framework is required. Download at http://gramian.de/emgr.m'); return; end

%%%%%%%% Setup %%%%%%%%

 N = 5;
 t = [0 0.01 2];
 q = [0 4*N 2*N];
 R = 4;
 p = ones(N,1);
 T = (t(3)-t(1))/t(2);

 X = [1.449;0.0; 0.400;-0.345; -1.125;0.448; -1.125;-0.448; 0.400;0.345; 0.0;-0.922; -1.335;0.810; -0.919;-0.349; 0.919;-0.349; 1.335;0.810]; %5-body eight-figure
 F = @(x,u,p) [x((2*N)+1:end);acc(x(1:2*N),u,p)];
 G = @(x,u,p)  x(1:2*N);

%%%%%%%% Original %%%%%%%%

 tic;
 Y = leapfrog(F,G,t,X,p);
 FULL = toc

%%%%%%%% Reduction %%%%%%%%

 tic;
 WO = emgr(F,G,q,t,'o',p,[0 0 0 0 0 0 0 0 0 0 0 2]);
 WOP = WO(1:(2*N),1:(2*N));
 WOV = WO((2*N)+1:end,(2*N)+1:end);
 [PP DD QQ] = svd(WOP); PP = PP(:,1:2*R); QQ = PP'; %diag(DD)'
 [TT DD VV] = svd(WOV); TT = TT(:,1:2*R); VV = TT'; %diag(DD)'
 OFFLINE = toc

%%%%%%%% Reduced %%%%%%%%

 f = @(x,u,p) [x((2*R)+1:end);QQ*acc(PP*x(1:2*R),u,p)];
 g = @(x,u,p) PP*x(1:2*R);
 x = [QQ*X(1:2*N);QQ*X((2*N)+1:end)];

 tic;
 y = leapfrog(f,g,t,x,p);
 ONLINE = toc

 ERROR = norm(norm(Y - y)./norm(Y))

%{
 f = @(x,u,p) [x((2*R)+1:end);VV*acc(TT*x(1:2*R),u,p)];
 g = @(x,u,p) TT*x(1:2*R);
 x = [VV*X(1:2*N);VV*X((2*N)+1:end)];

 tic;
 y = leapfrog(f,g,t,x,p);
 ONLINE = toc

 ERROR = norm(norm(Y - y)./norm(Y))
%}

%%%%%%%% Plot %%%%%%%%

 figure('PaperSize',[1.5,2.8],'PaperPosition',[0,0,2.8,1.5]);

 hold on;

 cmap = hsv(N+1);
 d = 1;

 for c=1:N
	plot(y(d,:),y(d+1,:),'--','Color',cmap(c,:));
	plot(y(d,end),y(d+1,end),'*','Color',cmap(c,:));
	d = d + 2;
 end

 ylim([-0.6 0.6]);
 hold off;

 if(nargin>0), print -dpng nbody.png; end

%%%%%%%% Acceleration %%%%%%%%

function y = acc(x,u,p)

 N = numel(x)/2;
 A = reshape(x,[2 N]);
 B = zeros(size(A));
 Z = zeros(1,N);
 y = zeros(2,N);

 for I=1:N
	B = bsxfun(@minus,A,A(:,I));
	Z = p'./(sqrt(0.0001+(sum(B.^2))).^3);
	B = bsxfun(@times,B,Z);
	y(:,I) = sum(B,2);
 end

 y = y(:);

%%%%%%%% Leapfrog %%%%%%%%

function y = leapfrog(f,g,t,z,p)

 T = (t(3)-t(1))/t(2);
 h = t(2);

 N = size(z,1);
 m = N/2;
 n = m + 1;
 l = f(z,0,p);
 k = l(n:N);

 for t=1:T
	z(1:m) = z(1:m) + h*z(n:N) + 0.5*h*h*k;
	l = f(z,0,p);
	z(n:N) = z(n:N) + 0.5*h*(k+l(n:N));
	y(:,t) = g(z,0,p);
	k = l(n:N);
 end
