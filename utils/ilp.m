function [A,B,C] = ilp(M,N,Q,s,r)
%%% summary: ilp (inverse lyapunov procedure)
%%% project: emgr - Empirical Gramian Framework ( http://gramian.de )
%%% authors: Christian Himpe ( 0000-0003-2194-6754 )
%%% license: 2-Clause BSD (2016--2017)
%$
    if(exist('emgr')~=2)
        disp('emgr framework is required. Download at http://gramian.de/emgr.m');
        return;
    end;

    if(nargin==3 || isempty(s))
        s = 0;
    end;

    if(nargin==5)
        rand('seed',r);
        randn('seed',r);
    end;

%% Gramian Eigenvalues
    WC = exp( 0.5*rand(N,1) );
    WO = exp( 0.5*rand(N,1) );

%% Gramian Eigenvectors
    [P,S,O] = svd(randn(N,N));

%% Balancing Transformation
    WC = P*diag(sqrt(WC))*P';
    WO = O*diag(sqrt(WO))*O';
    [U,D,V] = svd(WC*WO);

%% Input and Output
    B = randn(N,M);

    if(s)
        C = B';
    else
        C = randn(Q,N);
    end;

%% Scale Output Matrix
    BB = sum(B.*B,2);  % = diag(B*B')
    CC = sum(C.*C,1)'; % = diag(C'*C)
    SC = sqrt(BB./CC)';
    C = bsxfun(@times,C,SC);

%% Solve for System Matrix
    f = @(x,u,p,t) -D*x+B*u;
    g = @(x,u,p,t)  C*x;
    A = -emgr(f,g,[M,N,Q],[1.0/N,1.0],'c') - sqrt(eps)*speye(N);

%% Unbalance System
    if(s==0)
        A = V*A*U';
        B = V*B;
        C = C*U';
    end;
end
