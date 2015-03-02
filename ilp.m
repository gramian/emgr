function [A,B,C] = ilp(J,N,O,s,r)
% ilp (inverse lyapunov procedure)
% by Christian Himpe, 2013-2015 ( http://gramian.de )
% released under BSD 2-Clause License ( opensource.org/licenses/BSD-2-Clause )
%*

if(exist('emgr')~=2)
    disp('emgr framework is required. Download at http://gramian.de/emgr.m');
    return;
end

if(nargin==5)
    rand('seed',r);
    randn('seed',r);
end;

% Gramian Eigenvalues
 WC = exp( rand(N,1) );
 WO = exp( rand(N,1) );

% Gramian Eigenvectors
 [P S Q] = svd(randn(N));

% Balancing Transformation
 WC = P*diag(WC)*P';
 WO = Q*diag(WO)*Q';
 [U D V] = svd(WC*WO);

% Input and Output
 B = randn(N,J);

 if(nargin>=4 && s~=0),
        C = B';
 else,
        C = randn(O,N);
 end

% Scale Output Matrix
 BB = sum(B.*B,2);  % = diag(B*B')
 CC = sum(C.*C,1)'; % = diag(C'*C)
 SC = sqrt(BB./CC)';
 C = bsxfun(@times,C,SC);

% Solve System Matrix
 f = @(x,u,p) -0.01*D*x+B*u;
 g = @(x,u,p) C*x;
 A = -emgr(f,g,[J,N,O],[0,0.01,1.0],'c');

% Unbalance System
 A = V*A*U';
 B = V*B;
 C = C*U';

