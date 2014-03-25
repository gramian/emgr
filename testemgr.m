function testemgr(J)
% testemgr
% by Christian Himpe, 2013-2014 ( http://gramian.de )
% released under BSD 2-Clause License ( opensource.org/licenses/BSD-2-Clause )
%*

if(exist('emgr')~=2) disp('emgr framework is required. Download at http://gramian.de/emgr.m'); return; end

if(nargin<1) J=4; end	% number of inputs
N = J*J;		% number of states
O = J;			% number of outputs

S = 0;			% start time
h = 0.01;		% time step
T = 1.0;		% end time

A = rand(N,N);		% random system matrix
A(1:N+1:end) = -N;	% ensure stability
A = 0.5*(A+A');		% symmetrize system matrix
B = rand(N,J);		% random input matrix
C = B';			% ensure state-space symmetric system
P = A(:);		% parameter vector

f = @(x,u,p) reshape(p,[N N])*x+B*u;	% parametrized linear dynamic system
F = @(x,u,p) reshape(p,[N N])'*x+C'*u;	% adjoint dynamic system
g = @(x,u,p) C*x;			% linear output function 

G = -0.5*trace(C*inv(A)*B);		% System Gain	

disp('Computing Empirical Controllability Gramian (WC) and gain error: '); 
WC = emgr(f,g,[J N O],[S h T],'c',P);
abs(trace(WC)-G)

disp('Computing Empirical Observability Gramian (WO) and gain error: ');
WO = emgr(f,g,[J N O],[S h T],'o',P);
abs(trace(WO)-G)

disp('Computing Empirical Cross Gramian (WX) and gain error: ');
WX = emgr(f,g,[J N O],[S h T],'x',P);
abs(trace(WO)-G)

disp('Computing Empirical Approximate Cross Gramian (WY) and gain error: ');
WY = emgr(f,F,[J N O],[S h T],'y',P);
abs(trace(WY)-G)

disp('Computing Empirical Sensitivity Gramian (WS) and gain error: ');
WS = emgr(f,g,[J N O],[S h T],'s',P);
abs(trace(WS{1})-G)

disp('Computing Empirical Identifiability Gramian (WI) and gain error: ');
WI = emgr(f,g,[J N O],[S h T],'i',P);
abs(trace(WI{1})-G)

disp('Computing Empirical Joint Gramian (WJ) and gain error: ');
WJ = emgr(f,g,[J N O],[S h T],'j',P);
abs(trace(WJ{1})-G)

