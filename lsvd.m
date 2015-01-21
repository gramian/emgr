function [U,D,V] = lsvd(A,k,m)
% lsvd (lanczos sparse svd / truncated svd)
% by Christian Himpe, 2014-2015 ( http://gramian.de )
% released under BSD 2-Clause License ( opensource.org/licenses/BSD-2-Clause )
%
% Compatible with Matlab and Octave, for square matrices only!
%
% Based on:
%
% D. Bickson.
% "Gaussian Belief Propagation: Theory and Application".
% The Hebrew University of Jerusalem, 2008.
%
% J. Chen and Y. Saad.
% "Lanczos Vectors versus Singular Vectors for Effective Dimension Reduction".
% IEEE Transactions on Knowledge and Data Engineering, 21(8): 1091--1103, 2009.
%
%*
    if(nargin<3), m = k+1; end;

    n = size(A,2);

    a = zeros(m,1);
    b = zeros(m+1,1);
    q = zeros(n,1);
    r = zeros(n,1);

    Q = zeros(n,m+1);
    Q(:,1) = 0.5*ones(n,1);
    Q(:,1) = Q(:,1)/norm(Q(:,1));

    for I=1:m,

        q = Q(:,I);
        w = (A*(q'*A)') - b(I)*r;
        a(I) = w'*q;
        w = w - a(I)*q;

        w = w - Q(:,1:I-1)*(w'*Q(:,1:I-1))'; % re-orthogonalize

        b(I+1) = norm(w);
        Q(:,I+1) = w/b(I+1);
        r = q;
    end;

    T = diag(a);
    T(2:I+1:end) = b(2:end-1);
    T(I+1:I+1:end) = b(2:end-1);
    [u,d] = eig(T);
    [D,X] = sort(diag(d),'descend');

    [U,R] = qr(Q(:,1:I)*u(:,X(1:k)),0); % re-orthogonalize II
    D = sqrt(D(1:k));

    if(nargout==3)
        V = diag(1.0./D)*(U'*A);
    end
end

