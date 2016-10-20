function benchmark_non(o)
%%% summary: benchmark_non (nonlinear rc ladder benchmark)
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
    M = 1;
    N = 64;
    Q = M;
    T = [0.01,1.0];
    L = floor(T(2)/T(1)) + 1;
    U = @(t) ones(M,1)*(t<=0.5*L);
    X = zeros(N,1);

    g = @(x) exp(x)+x-1.0;

    A1 = spdiags(ones(N-1,1),-1,N,N)-speye(N);
    A2 = spdiags([ones(N-1,1);0],0,N,N)-spdiags(ones(N,1),1,N,N);

    NON = @(x,u,p,t) g(A1*x)-g(A2*x) + [u;sparse(N-1,1)];
    OUT = @(x,u,p,t) x(1);

%% FULL ORDER
    Y = ODE(NON,OUT,T,X,U,0);
    n1 = norm(Y(:),1);
    n2 = norm(Y(:),2);
    n8 = norm(Y(:),Inf);

%% OFFLINE
    tic;
    WX = emgr(NON,OUT,[M,N,Q],T,'x',0,0,U);
    [UU,D,VV] = svd(WX);
    OFFLINE = toc

%% EVALUATION
    for n=1:N-1
        uu = UU(:,1:n);
        vv = uu';
        x = vv*X;
        non = @(x,u,p,t) vv*NON(uu*x,u,p);
        out = @(x,u,p,t) OUT(uu*x,u,p);
        y = ODE(non,out,T,x,U,0);
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

