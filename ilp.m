function [A B C] = ilp(J,N,O,s)
% ilp (inverse lyapunov procedure)
% by Christian Himpe, 2013 ( http://gramian.de )
% released under BSD 2-Clause License ( http://gramian.de/#license )
%*

if(exist('emgr')~=2) disp('emgr framework is required. Download at http://gramian.de/emgr.m'); return; end

%% Gramian Eigenvalues
 WC = exp(-N + N*rand(N,1));
 WO = exp(-N + N*rand(N,1));

%% Gramian Eigenvectors
 X = randn(N,N);
 [U E V] = svd(X);

%% Balancing Trafo
 [P D Q] = svd(diag(WC.*WO));
 W = -D;

%% Input and Output
 B = randn(N,J);

 if(nargin<4 || s==0)
        C = randn(O,N);
 else
        C = B';
 end

%% Scale Output Matrix
 BB = sum(B.*B,2);  % = diag(B*B')
 CC = sum(C.*C,1)'; % = diag(C'*C)
 C = bsxfun(@times,C,sqrt(BB./CC)');

%% Solve System Matrix
 f = @(x,u,p) W*x+B*u;
 g = @(x,u,p) C*x;
 A = -emgr(f,g,[J N O],0,[0 0.01 1],'c');

%% Unbalance System
 T = U'*P';
 A = T*A*T';
 B = T*B;
 C = C*T';

