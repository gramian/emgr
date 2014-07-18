function [U,D] = lanczos(A,k)
% lanczos (lanczos sparse pod)
% by Christian Himpe, 2014 ( http://gramian.de )
% released under BSD 2-Clause License ( opensource.org/licenses/BSD-2-Clause )
%
% Compatible with Matlab and Octave, for square matrices only!
%
% Based on:
%
% D. Bickson, "Gaussian Belief Propagation: Theory and Application",
% The Hebrew University of Jerusalem, 2008
%
% J. Chen, Y. Saad, "Lanczos Vectors versus Singular Vectors for Effective Dimension Reduction",
% IEEE Transactions on Knowledge and Data Engineering, vol. 21, no. 8, pp. 1091-1103, August, 2009.
%
%*

n = size(A,2);
m = 2*k;

a = zeros(m,1); 
b = zeros(m+1,1); 
r = zeros(n,1);

Q = zeros(n,m+1); 
Q(:,1) = 0.5*rand(n,1);
Q(:,1) = Q(:,1)*(1.0/norm(Q(:,1),2));

p = -1;

for I=1:m

    q = Q(:,I);
    w = (A*(q'*A)') - b(I)*r;
    a(I) = w'*q;
    w = w - a(I)*q;

    w = w - Q(:,1:I-1)*(w'*Q(:,1:I-1))'; % re-orthogonalize

    b(I+1) = norm(w,2);
    Q(:,I+1) = w/b(I+1);
    r = q;

    if(I==m || (I>m && mod(I,20)==0) ), % starting at 2*k iterations, then every 20 test for convergence

        T = diag(a);
        T(2:I+1:end) = b(2:end-1);
        T(I+1:I+1:end) = b(2:end-1);
        [u,d] = eig(T);
        [D,X] = sort(diag(d),'descend');

        o = abs(sum(D(1:k)));
        if( abs(o-p) < 1e-8 || I>=n), break; end; %if sum of current and previous eigenvals is below tolerance
        p = o;

        m = m + 20;
        a(m,1) = 0;
        b(m+1,1) = 0;
        Q(n,m+1) = 0;
    end;

end;

U = Q(:,1:I)*u(:,X(1:k));
D = sqrt(D(1:k));

end
