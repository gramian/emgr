function benchmark_fss(o)
%%% summary: benchmark_fss (Flexible Space Structures Benchmark)
%%% project: emgr - Empirical Gramian Framework ( http://gramian.de )
%%% authors: Christian Himpe ( 0000-0003-2194-6754 )
%%% license: 2-Clause BSD (2017)
%$
    if(exist('emgr')~=2)
        error('emgr not found! Get emgr at: http://gramian.de');
    else
        global ODE;
        ODE = [];
        fprintf('emgr (version: %1.1f)\n',emgr('version'));
    end

    rand('seed',1009);

%% SETUPs

    K = 32;				% number of modes
    xi = rand(1,K)*0.001;		% damping ratio
    omega = rand(1,K)*100;		% natural frequencies

    N = 2*K;				% number of states
    M = 1;				% number of inputs
    Q = M;				% number of outputs

    A_k = cellfun(@(p) sparse([-2.0*p(1)*p(2),-p(2);p(2),0]), ...
                  num2cell([xi;omega],1),'UniformOutput',0);
    A = blkdiag(A_k{:});		% system matrix
    B = kron(rand(K,M),[1;0]);		% input matrix
    C = 10.0*rand(Q,2*K);		% output matrix

    h = 0.001;				% time step size
    T = 1.0;				% time horizon			
    L = floor(T/h) + 1;			% number of time steps
    U = @(t) ones(M,1)*(t<=h)/h;	% impulse input function
    X = zeros(N,1);

    LIN = @(x,u,p,t) A*x + B*u;		% vector field
    OUT = @(x,u,p,t) C*x;		% output functional

%% FULL ORDER
    Y = ODE(LIN,OUT,[h,T],X,U,0);
    %figure; plot(0:h:T,Y); return;
    n1 = norm(Y(:),1);
    n2 = norm(Y(:),2);
    n8 = norm(Y(:),Inf);

%% OFFLINE
    tic;
    s = 16;
    P = ceil(N/s);
    wx = cell(1,Q);
    for p=1:P
        wx{p} = emgr(LIN,OUT,[M,N,Q],[h,T],'x',0,[0,0,0,0,0,0,0,0,0,0,s,p]);
    end;
    [UU,DD,VV] = svd(cell2mat(wx));
    OFFLINE_TIME = toc

%% EVALUATION
    l1 = zeros(1,N-1);
    l2 = zeros(1,N-1);
    l8 = zeros(1,N-1);

    for n=1:N-1
        uu = UU(:,1:n);
        a = uu'*A*uu;
        b = uu'*B;
        c = C*uu;
        x = uu'*X;
        lin = @(x,u,p,t) a*x + b*u;
        out = @(x,u,p,t) c*x;
        y = ODE(lin,out,[h,T],x,U,0);
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

