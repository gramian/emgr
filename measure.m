function measure(o)
% nonlinearity measure
% by Christian Himpe, 2013 ( http://gramian.de )
% released under BSD 2-Clause License ( opensource.org/licenses/BSD-2-Clause )
%*

if(exist('emgr')~=2) disp('emgr framework is required. Download at http://gramian.de/emgr.m'); return; end

%%%%%%%% Setup %%%%%%%%

 J = 1;
 O = J;
 N = 16;
 T = [0 0.01 1];
 L = (T(3)-T(1))/T(2);
 X = zeros(N,1);

 rand('seed',1009);
 A = rand(N,N); A(1:N+1:end) = -0.55*N; A = 0.5*(A+A');
 B = rand(N,J);
 C = B';

 LIN = @(x,u,p) A*x + B*u;
 OUT = @(x,u,p) C*x;
 NIN = @(x,u,p) A*x + B*asinh(p*u);
 NST = @(x,u,p) A*asinh(p*x) + B*u;
 NOU = @(x,u,p) C*asinh(p*x);

%%%%%%%% Reduction %%%%%%%%%

% LINEAR
 WL = emgr(LIN,OUT,[J N O],0,T,'x');

% NONLINEAR
 K = 100;
 y = zeros(3,K);
 Q = 2.0/K;
 P = 0;

 for(I=1:K)
  Wi = emgr(NIN,OUT,[J N O],P,T,'c');
  Ws = emgr(NST,OUT,[J N O],P,T,'x');
  Wo = emgr(LIN,NOU,[J N O],P,T,'o');

  Ni = sum(sum(abs(WL-Wi)));
  Ns = sum(sum(abs(WL-Ws)));
  No = sum(sum(abs(WL-Wo)));

  y(:,I) = [Ni;Ns;No];
  P = P + Q;
 end

% NORMALIZE
 y = y./trace(WL);

%%%%%%%% Output %%%%%%%%

% TERMINAL
 I_S_O = sqrt(sum(abs(y).^2,2))

% PLOT
 if(nargin<1 || o==0 ), return; end
 l = (1:-0.01:0)'; cmap = [l,l,ones(101,1)];
 figure('PaperSize',[2.4,6.4],'PaperPosition',[0,0,6.4,2.4]);
 imagesc(y); caxis([0 max(max(y))]); colorbar; colormap(cmap); set(gca,'YTick',1:3); set(gca,'YTickLabel',{'I','S','O'}); set(gca,'XTickLabel',0:2/5:2);
 if(o==2 && exist('OCTAVE_VERSION')), print -dsvg measure.svg; end
