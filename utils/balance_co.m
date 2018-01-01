function [X,Y,Z] = balance_co(WC,WO)
%%% summary: balance_co (SVD-based balancing)
%%% project: emgr - EMpirical GRamian Framework ( http://gramian.de )
%%% authors: Christian Himpe ( 0000-0003-2194-6754 )
%%% license: 2-Clause BSD (2013--2017)
%$
    if(nargin==1)
        [VL,EL] = eig(WC);
        EB = VL' \ VR;
        [UB,DB,VB] = svd(EB);
        HSV = diag(DB);
        db = sqrt(1.0./DB));
        TL = (VL*UB*diag(db))';
        TR = VR*VB'*diag(db);
    else if(nargin==2)
        [L,D,l] = svd(WC);
        LC = L*diag(sqrt(diag(D)));
        [L,D,l] = svd(WO);
        LO = L*diag(sqrt(diag(D)));
        [U,Y,V] = svd(LO'*LC);
        X = ( LO*U*diag(1.0./sqrt(diag(Y))) )';
        Z =   LC*V*diag(1.0./sqrt(diag(Y)));
    end
end
