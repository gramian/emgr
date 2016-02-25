function state_bt(o)
% state_bt (balanced truncation linear state reduction)
% by Christian Himpe, 2013-2016 ( http://gramian.de )
% released under BSD 2-Clause License ( opensource.org/licenses/BSD-2-Clause )
%*
    if(exist('emgr')~=2)
        error('emgr not found! Get emgr at: http://gramian.de');
    else
        global ODE; ODE = [];
        fprintf('emgr (version: %g)\n',emgr('version'));
    end

%% SETUP
    J = 8;
    N = 64;
    O = J;
    T = [0.01,1.0];
    L = floor(T(2)/T(1)) + 1;
    U = [ones(J,1),zeros(J,L-1)];
    X = zeros(N,1);

    rand('seed',1009);
    A = rand(N,N);
    A(1:N+1:end) = -0.55*N;
    A = 0.5*(A+A');
    B = rand(N,J);
    C = B';

    LIN = @(x,u,p) A*x + B*u;
    OUT = @(x,u,p) C*x;

%% FULL ORDER
    Y = ODE(LIN,OUT,T,X,U,0);
    n1 = norm(Y(:),1);
    n2 = norm(Y(:),2);
    n8 = norm(Y(:),Inf);

%% OFFLINE
    tic;
    WC = emgr(LIN,OUT,[J,N,O],T,'c');
    WO = emgr(LIN,OUT,[J,N,O],T,'o');
    [VV,D,UU] = balance(WC,WO);
    OFFLINE = toc

%% EVALUATION
    for I=1:N-1
        uu = UU(:,1:I);
        vv = VV(1:I,:);
        a = vv*A*uu;
        b = vv*B;
        c = C*uu;
        x = vv*X;
        lin = @(x,u,p) a*x + b*u;
        out = @(x,u,p) c*x;
        y = ODE(lin,out,T,x,U,0);
        l1(I) = norm(Y(:)-y(:),1)/n1;
        l2(I) = norm(Y(:)-y(:),2)/n2;
        l8(I) = norm(Y(:)-y(:),Inf)/n8;
    end;

%% OUTPUT
    if(nargin>0 && o==0), return; end; 
    figure('Name',mfilename,'NumberTitle','off');
    semilogy(1:N-1,[l1;l2;l8],{'r','g','b'},'linewidth',2);
    xlim([1,N-1]);
    ylim([10^floor(log10(min([l1(:);l2(:);l8(:)]))-1),1]);
    pbaspect([2,1,1]);
    legend('L1 Error ','L2 Error ','L8 Error ','location','northeast');
    if(nargin>0 && o==1), print('-dsvg',[mfilename(),'.svg']); end;
end

%% ======== Balancer ========
function [X,Y,Z] = balance(WC,WO)

    [L,D,l] = svd(WC); LC = L*diag(sqrt(diag(D)));
    [L,D,l] = svd(WO); LO = L*diag(sqrt(diag(D)));
    [U,Y,V] = svd(LO'*LC);
    X = ( LO*U*diag(1.0./sqrt(diag(Y))) )';
    Z =   LC*V*diag(1.0./sqrt(diag(Y)));
end
