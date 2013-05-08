function measure(o)
% measure (nonlinearity measure)
% by Christian Himpe, 2013 ( http://gramian.de )
% released under BSD 2-Clause License ( http://gramian.de/#license )
%*

if(exist('emgr')~=2) disp('emgr framework is required. Download at http://gramian.de/emgr.m'); return; end

%%%%%%%% Setup %%%%%%%%

%% CONSTANTS
 J = 1;
 O = J;
 N = 16;
 R = O;
 t = [0 0.01 1];
 T = (t(3)-t(1))/t(2);
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
 LIN = @(x,u,p) A*x + B*u;
 OUT = @(x,u,p) C*x;
 NIN = @(x,u,p) A*x + B*asinh(p*u);
 NSY = @(x,u,p) A*asinh(p*x) + B*u;
 NOU = @(x,u,p) C*asinh(p*x);
%%

%%%%%%%% Linear State Reduction %%%%%%%%%

%% WX (linear)
 Wl = emgr(LIN,OUT,[J N O],0,t,'x');
%%

K = 100;
y = zeros(3,K);
q = 2/K;
p = q;

for(I=1:K)

%% WX (nonlinear)
 Wi = emgr(NIN,OUT,[J N O],p,t,'x');
 Ws = emgr(NSY,OUT,[J N O],p,t,'x');
 Wo = emgr(LIN,NOU,[J N O],p,t,'x');

 Ni = sum(sum(abs(Wl-Wi)));
 Ns = sum(sum(abs(Wl-Ws)));
 No = sum(sum(abs(Wl-Wo)));

 y(:,I) = [Ni;Ns;No];
 p = p + q;
%%

end

%% Normalize
Nl = sum(sum(Wl));
y = y./Nl;
I_S_O = sqrt(sum(abs(y).^2,2))
%%

%% PLOT
 if(nargin<1 || o==0 ) return; end
 l = (1:-0.01:0)'; cmap = [l,l,ones(101,1)];
 figure('PaperSize',[2.4,6.4],'PaperPosition',[0,0,6.4,2.4]);
 imagesc(y); caxis([0 max(max(y))]); colorbar; colormap(cmap);
 set(gca,'YTick',1:3);
 set(gca,'YTickLabel',{'I','S','O'});
 set(gca,'XTickLabel',0:2/5:2);
 if(o==2 &&  exist('OCTAVE_VERSION'))
 	print -dsvg measure.svg;
 end
%%
