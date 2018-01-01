function benchmark_lin(o)
%%% summary: benchmark_lin (ISS 1-r model reduction benchmark)
%%% project: emgr - EMpirical GRamian Framework ( http://gramian.de )
%%% authors: Christian Himpe ( 0000-0003-2194-6754 )
%%% license: 2-Clause BSD (2013--2017)
%$
    if(exist('emgr')~=2)
        error('emgr not found! Get emgr at: http://gramian.de');
    else
        global ODE;
        ODE = [];
        fprintf('emgr (version: %1.1f)\n',emgr('version'));
    end

    D = 'iss';
    try
        if(exist(['/tmp/',D,'.mat'],'file')==0)
            urlwrite(['http://slicot.org/objects/software/shared/bench-data/',D,'.zip'],['/tmp/',D,'.zip']);
        end
    catch
        disp('Cannot download benchmark data!');
        return;
    end
    unzip(['/tmp/',D,'.zip'],'/tmp');
    load(['/tmp/',D,'.mat']);

%% SYSTEM SETUP
    N = size(A,2);			% number of states
    K = N/2;				% number of coupled states
    M = size(B,2);			% number of inputs
    Q = size(C,1);			% number of outputs
    h = 0.01;				% time step size
    T = 1.0;				% time horizon
    X = zeros(N,1);			% initial state
    U = @(t) ones(M,1)*(t<=h)/h;	% impulse input function

    LIN = @(x,u,p,t) A*x + B*u;		% vector field
    OUT = @(x,u,p,t) C*x;		% output functional

%% FULL ORDER MODEL REFERENCE SOLUTION
    Y = ODE(LIN,OUT,[h,T],X,U,0);
    n1 = norm(Y(:),1);
    n2 = norm(Y(:),2);
    n8 = norm(Y(:),Inf);

%% STRUCTURED REDUCED ORDER MODEL PROJECTION ASSEMBLY
    tic;
    WX = emgr(LIN,OUT,[M,N,Q],[h,T],'x'); % use velocity cross gramian
    WXV = WX(K+1:N,K+1:N);
    [UU,D,VV] = svd(WXV);
    OFFLINE_TIME = toc

    %{
    tic;
    WX = emgr(LIN,OUT,[M,N,Q],[h,T],'x'); % use position cross gramian
    WXP = WX(1:K,1:K);
    [UU,D,VV] = svd(WXP);
    OFFLINE_TIME = toc
    %}

%% REDUCED ORDER MODEL EVALUATION
    l1 = zeros(1,K-1);
    l2 = zeros(1,K-1);
    l8 = zeros(1,K-1);

    for n=1:K-1
        uu = UU(:,1:n);
        vv = uu';
        a = [zeros(n,n),eye(n);vv*A(K+1:N,1:K)*uu,vv*A(K+1:N,K+1:N)*uu];
        b = [zeros(n,M);vv*B(K+1:N,:)];
        c = [zeros(Q,n),C(:,K+1:N)*uu];
        x = zeros(2*n,1);
        lin = @(x,u,p,t) a*x + b*u;
        out = @(x,u,p,t) c*x;
        y = ODE(lin,out,[h,T],x,U,0);
        l1(n) = norm(Y(:)-y(:),1)/n1;
        l2(n) = norm(Y(:)-y(:),2)/n2;
        l8(n) = norm(Y(:)-y(:),Inf)/n8;
    end;

%% PLOT REDUCED ORDER VS RELATIVE ERRORS
    if(nargin>0 && o==0), return; end; 
    figure('Name',mfilename,'NumberTitle','off');
    semilogy(2:2:N-2,l1,'r','linewidth',2); hold on;
    semilogy(2:2:N-2,l2,'g','linewidth',2);
    semilogy(2:2:N-2,l8,'b','linewidth',2); hold off;
    xlim([2,N-2]);
    ylim([1e-16,1]);
    pbaspect([2,1,1]);
    legend('L1 Error ','L2 Error ','L8 Error ','location','northeast');
    set(gca,'YGrid','on');
    if(nargin>0 && o==1), print('-dsvg',[mfilename(),'.svg']); end;
end

