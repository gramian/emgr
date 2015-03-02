function vernval(J)
% vernval (verification and validation)
% by Christian Himpe, 2013-2015 ( http://gramian.de )
% released under BSD 2-Clause License ( opensource.org/licenses/BSD-2-Clause )
%*

if(exist('emgr')~=2)
    disp('emgr framework is required. Download at http://gramian.de/emgr.m');
    return;
end

rand('seed',1009);

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
P = ones(N,1);		% parameter vector

f = @(x,u,p) A*x+B*u;   % linear dynamic system vector field
g = @(x,u,p) C*x;       % linear output functional

G = @(x,u,p) A'*x+C'*u;	% adjoint dynamic system vector field

F = @(x,u,p) A*x+B*u+p; % linear parametrized system vector field


WC = emgr(f,g,[J,N,O],[S,h,T],'c');
fprintf('Empirical Controllability Gramian: WC - OK\n');

Wc = emgr(G,g,[J,N,O],[S,h,T],'c');
fprintf('Adjoint Empirical Controllability Gramian: Wc - OK\n');

WO = emgr(f,g,[J,N,O],[S,h,T],'o');
fprintf('Empirical Observability Gramian: WO - OK\n');

WX = emgr(f,g,[J,N,O],[S,h,T],'x');
fprintf('Empirical Cross Gramian: WX - OK\n');

WY = emgr(f,G,[J,N,O],[S,h,T],'y');
fprintf('Empirical Linear Cross Gramian: WY - OK\n');

dWCWO = norm(WC-WO,'fro')
dWcWO = norm(Wc-WO,'fro')
dWCWY = norm(WC-WY,'fro')
dWCWX = norm(WC-WX,'fro')
dWOWX = norm(WO-WX,'fro')
dWXWY = norm(WX-WY,'fro')


WC = emgr(F,g,[J,N,O],[S,h,T],'c',P);
fprintf('Parametrized Empirical Controllability Gramian: WC - OK\n');

WO = emgr(F,g,[J,N,O],[S,h,T],'o',P);
fprintf('Parametrized Empirical Observability Gramian: WO - OK\n');

WX = emgr(F,g,[J,N,O],[S,h,T],'x',P);
fprintf('Parametrized Empirical Cross Gramian: WX - OK\n');

WS = emgr(F,g,[J,N,O],[S,h,T],'s',P);
fprintf('Empirical Sensitivity Gramian: WS - OK\n');

WI = emgr(F,g,[J,N,O],[S,h,T],'i',P);
fprintf('Empirical Identifiability Gramian: WI - OK\n');

WJ = emgr(F,g,[J,N,O],[S,h,T],'j',P);
fprintf('Empirical Joint Gramian: WJ - OK\n');

dWCWS = norm(WC-WS{1},'fro')
dWOWI = norm(WO-WI{1},'fro')
dWXWJ = norm(WX-WJ{1},'fro')
dWIWJ = norm(WI{2}-WJ{2},'fro')

