function [X,Y,Z] = balance_co(WC,WO,o)
% balance_co (balance controllability and observability gramians)
% by Christian Himpe, 2014-2015 ( http://gramian.de )
% released under BSD 2-Clause License ( opensource.org/licenses/BSD-2-Clause )
%*
    if(nargin<1 || o==0)
        [L,D,l] = svd(WC);
        LC = L*diag(sqrt(diag(D)));

        [L,D,l] = svd(WO);
        LO = L*diag(sqrt(diag(D)));
    else
        LC = chol(WC);

        LO = chol(WO);
    end

    [U,Y,V] = svd(LO'*LC);

    % Y : Hankel Singular Values
    X = ( LO*U*diag(1.0./sqrt(diag(Y))) )';	% Left Projection
    Z =   LC*V*diag(1.0./sqrt(diag(Y)));	% Right Projection
end
