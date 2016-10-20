function test_dwjz(o)
%%% summary: test_dwjz (distributed cross-identifiability gramian reduction)
%%% project: emgr - Empirical Gramian Framework ( http://gramian.de )
%%% authors: Christian Himpe ( 0000-0003-2194-6754 )
%%% license: 2-Clause BSD (2016)
%$
    if(exist('emgr')~=2)
        error('emgr not found! Get emgr at: http://gramian.de');
    else
        global ODE;
        ODE = [];
        fprintf('emgr (version: %1.1f)\n',emgr('version'));
    end

%% SETUP
    M = 4;
    N = M*M*M;
    Q = M;
    T = [0.01,1.0];
    L = floor(T(2)/T(1)) + 1;
    U = @(t) ones(M,1)*(t<=T(1))/T(1);
    X = zeros(N,1);
    P = 0.5+0.5*cos(1:N)';

    A = -gallery('lehmer',N);
    B = toeplitz(1:N,1:M)./N;
    C = B';

    LIN = @(x,u,p,t) A*x + B*u + p;
    OUT = @(x,u,p,t) C*x;

%% FULL ORDER
    Y = ODE(LIN,OUT,T,X,U,P);
    %figure; plot(0:T(1):T(2),Y); return;
    n1 = norm(Y(:),1);
    n2 = norm(Y(:),2);
    n8 = norm(Y(:),Inf);

%% OFFLINE
    global DWX;
    DWX = [];

    tic;
    WJ = emgr(LIN,OUT,[M,N,Q],T,'j',[zeros(N,1),ones(N,1)],[0,0,0,0,0,0,1,0,0,0]);
    [UU,D,VV] = svd(WJ{2});
    OFFLINE = toc

    tic;

    w = 8;
    K = ceil(2*N/w);
    wj = [];
    for k=1:K
        DWX = [w,k];
        wj = [wj,emgr(LIN,OUT,[M,N,Q],T,'j',[zeros(N,1),ones(N,1)],[0,0,0,0,0,0,1,0,0,0])];
    end;
    wx = wj(:,1:N);
    wi = -0.5*wj(:,N+1:end)'*ainv(wx+wx')*wj(:,N+1:end);

    OFFLINE = toc
    RESIDUAL1 = norm(WJ{1}-wx)
    RESIDUAL2 = norm(WJ{2}-wi)
    DWX = [];

%% EVALUATION
    for n=1:N-1
        uu = UU(:,1:n);
        y = ODE(LIN,OUT,T,X,U,uu*uu'*P);
        l1(n) = norm(Y(:)-y(:),1)/n1;
        l2(n) = norm(Y(:)-y(:),2)/n2;
        l8(n) = norm(Y(:)-y(:),Inf)/n8;
    end;

%% OUTPUT
    if(nargin>0 && o==0), return; end; 
    figure('Name',mfilename,'NumberTitle','off');
    semilogy(1:N-1,l1,'r','linewidth',2); hold on;
    semilogy(1:N-1,l2,'g','linewidth',2);
    semilogy(1:N-1,l8,'b','linewidth',2); hold off;
    xlim([1,N-1]);
    ylim([1e-16,1]);
    pbaspect([2,1,1]);
    legend('L1 Error ','L2 Error ','L8 Error ','location','southeast');
    if(nargin>0 && o==1), print('-dsvg',[mfilename(),'.svg']); end;
end

function x = ainv(m) 
%%% summary: ainv (approximate inverse)
%$
    d = diag(m);
    d(d~=0) = 1.0./d(d~=0);
    n = numel(d);
    x = bsxfun(@times,m,-d);
    x = bsxfun(@times,x,d');
    x(1:n+1:end) = d;
end

