function [U,D,V] = lanczos(A,k,m)
% lanczos (lanczos sparse svd)
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

if(nargin==2), m = 2*k;

n = size(A,2);
m = min(max(m,k),n);

a = zeros(m,1);
b = zeros(m+1,1);
r = zeros(n,1);

Q = zeros(n,m+1);
q = 0.5*ones(n,1);
q = q*(1.0/norm(q(:,1),2));
Q(:,1) = q;

for I=1:m

    q = Q(:,I);
    w = (A*(q'*A)') - b(I)*r;
    a(I) = w'*q;
    w = w - a(I)*q;

    w = w - Q(:,1:I-1)*(w'*Q(:,1:I-1))';

    b(I+1) = norm(w,2);
    Q(:,I+1) = w*(1.0/b(I+1));
    r = q;

end;

T = diag(a); 
T(2:I+1:end) = b(2:end-1);
T(I+1:I+1:end) = b(2:end-1);
[u,d] = eig(T);
[D,X] = sort(diag(d),'descend');

U = Q(:,1:I)*u(:,X(1:k));
D = sqrt(D(1:k));

if(nargout==3),
    V = zeros(k,n);
    for I=1:k
        V(I,:) = (U(:,I)'*A)/D(I);
    end;
end;

end
