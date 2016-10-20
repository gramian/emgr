function benchmark_lin(o)
%%% summary: benchmark_lin (iss model reduction benchmark)
%%% project: emgr - Empirical Gramian Framework ( http://gramian.de )
%%% authors: Christian Himpe ( 0000-0003-2194-6754 )
%%% license: 2-Clause BSD (2013--2016)
%$
    if(exist('emgr')~=2)
        error('emgr not found! Get emgr at: http://gramian.de');
    else
        global ODE;
        ODE = [];
        fprintf('emgr (version: %1.1f)\n',emgr('version'));
    end

    D = 'iss';
    if(exist(['/tmp/',D,'.mat'],'file')==0)
        urlwrite(['http://slicot.org/objects/software/shared/bench-data/',D,'.zip'],['/tmp/',D,'.zip']);
        unzip(['/tmp/',D,'.zip'],'/tmp');
    end
    load(['/tmp/',D,'.mat']);

%% SETUP
    N = size(A,2);
    K = N/2;
    M = size(B,2);
    Q = size(C,1);
    T = [0.01,1.0];
    L = floor(T(2)/T(1)) + 1;
    U = @(t) ones(M,1)*(t<=T(1))/T(1);
    X = zeros(N,1);

    LIN = @(x,u,p,t) A*x + B*u;
    OUT = @(x,u,p,t) C*x;

%% FULL ORDER
    Y = ODE(LIN,OUT,T,X,U,0);
    n1 = norm(Y(:),1);
    n2 = norm(Y(:),2);
    n8 = norm(Y(:),Inf);

%% OFFLINE
    tic;
    WX = emgr(LIN,OUT,[M,N,Q],T,'x'); % use velocity cross gramian
    WXV = WX(K+1:N,K+1:N);
    [UU,D,VV] = svd(WXV);
    OFFLINE = toc

    %{
    tic;
    WX = emgr(LIN,OUT,[M,N,Q],T,'x'); % use position cross gramian
    WXP = WX(1:K,1:K);
    [UU,D,VV] = svd(WXP);
    OFFLINE = toc
    %}

%% EVALUATION
    for n=1:K-1
        uu = UU(:,1:n);
        vv = uu';
        a = [zeros(n,n),eye(n);vv*A(K+1:N,1:K)*uu,vv*A(K+1:N,K+1:N)*uu];
        b = [zeros(n,M);vv*B(K+1:N,:)];
        c = [zeros(Q,n),C(:,K+1:N)*uu];
        x = zeros(2*n,1);
        lin = @(x,u,p,t) a*x + b*u;
        out = @(x,u,p,t) c*x;
        y = ODE(lin,out,T,x,U,0);
        l1(n) = norm(Y(:)-y(:),1)/n1;
        l2(n) = norm(Y(:)-y(:),2)/n2;
        l8(n) = norm(Y(:)-y(:),Inf)/n8;
    end;

%% OUTPUT
    if(nargin>0 && o==0), return; end; 
    figure('Name',mfilename,'NumberTitle','off');
    semilogy(2:2:N-2,l1,'r','linewidth',2); hold on;
    semilogy(2:2:N-2,l2,'g','linewidth',2);
    semilogy(2:2:N-2,l8,'b','linewidth',2); hold off;
    xlim([2,N-2]);
    ylim([10^floor(log10(min([l1(:);l2(:);l8(:)]))-1),1]);
    pbaspect([2,1,1]);
    legend('L1 Error ','L2 Error ','L8 Error ','location','northeast');
    if(nargin>0 && o==1), print('-dsvg',[mfilename(),'.svg']); end;
end

