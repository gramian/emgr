function testemgr(N,J)
% testemgr
% by Christian Himpe, 2013 ( http://gramian.de )
% released under BSD 2-Clause License ( opensource.org/licenses/BSD-2-Clause )
%*

if(exist('emgr')~=2) disp('emgr framework is required. Download at http://gramian.de/emgr.m'); return; end

if(nargin<1) N=8; end	%number of states
if(nargin<2) J=4; end	%number of inputs
O = J;			%number of outputs

S = 0;			%start time
h = 0.01;		%time step
T = 1;			%end time

A = rand(N,N);		%random system matrix
A(1:N+1:end) = -N;	%ensure stability
A = 0.5*(A+A');		%symmetrize system matrix
B = rand(N,J);		%random input matrix
C = B';			%ensure symmetric system

f = @(x,u,p) reshape(p,[N N])*x+B*u;	%parametrized linear dynamic system
g = @(x,u,p) C*x;			%linear output function 


disp('Computing Empirical Controllability Gramian (WC)');

WC = emgr(f,g,[J N O],A(:),[S h T],'c');

disp('Computing Empirical Observability Gramian (WO) ');

WO = emgr(f,g,[J N O],A(:),[S h T],'o');

disp('Computing Empirical Cross Gramian (WX)');

WX = emgr(f,g,[J N O],A(:),[S h T],'x');

disp('Computing Balanced Proper Orthogonal Decomposition (BPOD)');

F = @(x,u,p) reshape(p,[N N])*x+C'*u;	%dual system
G = @(x,u,p) B'*x;			%dual output

BC = emgr(f,g,[J N O],A(:),[S h T],'c');
BO = emgr(F,G,[O N J],A(:),[S h T],'c');

disp('Computing Empirical Covariance Matrices');

u = [1;2;3;4]*(1./[0.01:0.1:10]);		%custom input

VC = emgr(f,g,[J N O],A(:),[S h T],'c',0,u);
VO = emgr(f,g,[J N O],A(:),[S h T],'o',0,u);
VX = emgr(f,g,[J N O],A(:),[S h T],'x',0,u);


fprintf('\nSquareroot of Eigenvalues of WC*WO: ');	lambda_WCWO = sort(sqrt(eig(WC*WO)),'descend')'

fprintf('Absolute Value of Eigenvalues of WX: '); lambda_WX   = sort(abs(eig(WX)),'descend')'

fprintf('Squareroot of Eigenvalues of Balanced POD: '); lambda_BPOD = sort(sqrt(eig(BO*BC)),'descend')'

fprintf('Trace of WX:               '); Tr_WX = trace(WX)

fprintf('Trace of Half System Gain: '); Tr_Gn = trace(-0.5*C*inv(A)*B)

W = emgr(f,g,[J N O],A(:),[S h T],'s');
Wc = W{1};				%controllability gramian
WS = W{2};				%parameter controllability gramian
fprintf('Error in byproduct WC from WS: '); E_WCWc = max(max(abs(WC-Wc)))

W = emgr(f,g,[J N O],A(:),[S h T],'i');
Wo = W{1};				%observability gramian
WI = W{2};				%parameter observability gramian
fprintf('Error in byproduct WO from WI: '); E_WOWo = max(max(abs(WO-Wo)))

W = emgr(f,g,[J N O],A(:),[S h T],'j');
Wx = W{1};				%cross gramian
WJ = W{2};				%parameter cross identifiability gramian
fprintf('Error in byproduct WX from WJ: '); E_WXWx = max(max(abs(WX-Wx)))



