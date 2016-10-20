function combined_wj(o)
%%% summary: combined_wj (joint gramian nonlinear combined reduction)
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
    U = @(t) ones(M,1)*(t<=T(1))/T(1);
    X = zeros(N,1);

    rand('seed',1009);

    A = full(sprand(N,N,1./N));
    A(1:N+1:end) = -2.0;
    A = A * 10.0;
    B = 10.0*full(sprand(N,M,2.0/N));
    C = rand(Q,N);

    P = 0.9*rand(N,10) + 0.1;
    R = [0.1*ones(N,1),ones(N,1)];

    NON = @(x,u,p,t) A*tanh(p.*x) + B*u;
    OUT = @(x,u,p,t) C*x;

%% OFFLINE
    tic;
    WJ = emgr(NON,OUT,[M,N,Q],T,'j',R,[0,0,0,0,0,0,0,1,1,0]);
    [UU,D,VV] = svd(WJ{1});
    [PP,D,QQ] = svd(WJ{2});
    OFFLINE = toc

%% EVALUATION
    for H=1:1
        Y = ODE(NON,OUT,T,X,U,P(:,H));
        n2 = norm(Y(:),2);

        a = 1;
        for n=1:4:N
            uu = UU(:,1:n);
            vv = uu';
            b = 1;
            for K=1:4:N
                pp = PP(:,1:K);
                qq = pp';
                p = pp*qq*P(:,H);
                x = vv*X;
                non = @(x,u,p,t) vv*NON(uu*x,u,p);
                out = @(x,u,p,t) OUT(uu*x,u,p);
                y = ODE(non,out,T,x,U,p);
                l2(a,b,H) = min(1.0,norm(Y(:)-y(:),2)/n2);
                b = b + 1;
            end;
            a = a + 1;
        end;
    end;
    l2 = sqrt(sum(l2.^2,3));

%% OUTPUT
    if(nargin>0 && o==0), return; end; 
    figure('Name',mfilename,'NumberTitle','off');
    h = surf(l2);
    xlim([1,16]);
    ylim([1,16]);
    set(gca,'ZScale','log');
    set(gca,'XTick',1:2:16,'XTickLabel',1:8:N);
    set(gca,'YTick',1:2:16,'YTickLabel',1:8:N);
    ylabel('State Dimension')
    xlabel('Parameter Dimension');
    if(exist('viridis')==0), colormap(hot); end;
    set(h,'CData',log10(get(h,'CData')))
    view(135,30);
    if(nargin>0 && o==1), print('-dsvg',[mfilename(),'.svg']); end;
end

