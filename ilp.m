function [A,B,C] = ilp(J,N,O,s,r)
% ilp (inverse lyapunov procedure)
% by Christian Himpe, 2013-2015 ( http://gramian.de )
% released under BSD 2-Clause License ( opensource.org/licenses/BSD-2-Clause )
%*
    if(exist('emgr')~=2)
        error('emgr not found! Get emgr at: http://gramian.de');
    end

    if(nargin==5)
        rand('seed',r);
        randn('seed',r);
    end;

    % Gramian Eigenvalues
    wc = exp( rand(N,1) );
    wo = exp( rand(N,1) );

    % Gramian Eigenvectors
    [P,S,Q] = svd(randn(N));

    % Balancing Transformation
    WC = P*diag(wc)*P';
    WO = Q'*diag(wo)*Q;
    [U,D,V] = svd(WC*WO);

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
    f = @(x,u,p) -D*x+B*u;
    g = @(x,u,p) C*x;
    A = -emgr(f,g,[J,N,O],[1.0/N,1.0],'c') - (1.0/N)*speye(N);

    % Unbalance System
    A = V'*A*U;
    B = V'*B;
    C = C*U;
end
