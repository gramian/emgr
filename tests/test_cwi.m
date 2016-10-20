function test_cwi(o)
%%% summary: test_cwi (identifiability gramian combined reduction)
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
    tic;
    WC = emgr(LIN,OUT,[M,N,Q],T,'c');    
    WI = emgr(LIN,OUT,[M,N,Q],T,'i',[zeros(N,1),ones(N,1)]);
    [L1,D1,R1] = svd(WC); LC = L1*diag(sqrt(diag(D1)));
    [L2,D2,R2] = svd(WI{1}); LO = L2*diag(sqrt(diag(D2)));
    [Lb,Db,Rb] = svd(LO'*LC);
    UU = ( LO*Lb*diag(1.0./sqrt(diag(Db))) )';
    VV =   LC*Rb*diag(1.0./sqrt(diag(Db)));
    [PP,D,QQ] = svd(WI{2});
    OFFLINE = toc

%% EVALUATION
    for n=1:N-1
        uu = VV(:,1:n);
        vv = UU(1:n,:);
        pp = PP(:,1:n);
        lin = @(x,u,p,t) vv*LIN(uu*x,u,p);
        out = @(x,u,p,t) OUT(uu*x,u,p);
        y = ODE(lin,out,T,vv*X,U,pp*pp'*P);
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
    legend('L1 Error ','L2 Error ','L8 Error ','location','northeast');
    if(nargin>0 && o==1), print('-dsvg',[mfilename(),'.svg']); end;
end

