function energy_wz(o)
%%% summary: energy_wz (cross gramian experimental reduction)
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

%% SETUP
    M = 4;
    N = M*M*M;
    Q = 1;
    T = [0.01,1.0];
    L = floor(T(2)/T(1)) + 1;
    U = @(t) ones(M,1)*(t<=T(1))/T(1);
    X = zeros(N,1);

    rand('seed',1009);
    A = rand(N,N);
    A(1:N+1:end) = -0.55*N;
    B = rand(N,M);
    C = rand(Q,N);

    LIN = @(x,u,p,t) A*x + B*u;
    OUT = @(x,u,p,t) x'*x;

%% FULL ORDER
    Y = ODE(LIN,OUT,T,X,U,0);
    n1 = norm(Y(:),1);
    n2 = norm(Y(:),2);
    n8 = norm(Y(:),Inf);

% OFFLINE
    tic;
    WZ = emgr(LIN,OUT,[M,N,Q],T,'x',0,[0,0,0,0,0,0,1,0,0,0],1,0,1);
    [UU,D,VV] = svd(WZ);
    OFFLINE = toc

%% EVALUATION
    for n=1:N-1
        uu = UU(:,1:n);
        vv = uu';
        a = vv*A*uu;
        b = vv*B;
        x = vv*X;
        lin = @(x,u,p,t) a*x + b*u;
        out = @(x,u,p,t) OUT(uu*x,u,p);
        y = ODE(lin,out,T,x,U,0);
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
    ylim([10^floor(log10(min([l1(:);l2(:);l8(:)]))-1),1]);
    pbaspect([2,1,1]);
    legend('L1 Error ','L2 Error ','L8 Error ','location','northeast');
    if(nargin>0 && o==1), print('-dsvg',[mfilename(),'.svg']); end;
end

