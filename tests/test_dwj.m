function test_dwj(o)
%%% summary: test_dwj (distributed cross-identifiability gramian reduction)
%%% project: emgr - Empirical Gramian Framework ( http://gramian.de )
%%% authors: Christian Himpe ( 0000-0003-2194-6754 )
%%% license: 2-Clause BSD (2014--2017)
%$
    if(exist('emgr')~=2)
        error('emgr not found! Get emgr at: http://gramian.de');
    else
        global ODE;
        ODE = [];
        fprintf('emgr (version: %1.1f)\n',emgr('version'));
    end

%% SYSTEM SETUP
    M = 4;				% number of inputs
    N = M*M*M;				% number of states
    Q = M;				% number of outputs
    h = 0.01;				% time step size
    T = 1.0;				% time horizon
    X = zeros(N,1);			% initial state
    U = @(t) ones(M,1)*(t<=h)/h;	% impulse input function
    P = 0.5+0.5*cos(1:N)';		% parameter
    R = [zeros(N,1),ones(N,1)];		% parameter range

    A = -gallery('lehmer',N);		% system matrix
    B = toeplitz(1:N,1:M)./N;		% input matrix
    C = B';				% output matrix

    LIN = @(x,u,p,t) A*x + B*u + p;	% vector field
    OUT = @(x,u,p,t) C*x;		% output functional

%% FULL ORDER MODEL REFERENCE SOLUTION
    Y = ODE(LIN,OUT,[h,T],X,U,P);
    %figure; plot(0:h:T,Y); return;
    n1 = norm(Y(:),1);
    n2 = norm(Y(:),2);
    n8 = norm(Y(:),Inf);

%% COMPARATIVE REDUCED ORDER MODEL PROJECTION ASSEMBLY
    tic;
    WJ = emgr(LIN,OUT,[M,N,Q],[h,T],'j',R);
    [UU,D,VV] = svd(WJ{2});
    OFFLINE_TIME_FULL = toc

    tic;

    w = 8;
    K = ceil(2*N/w);
    wj = cell(1,K);
    for k=1:K
        wj{k} = emgr(LIN,OUT,[M,N,Q],[h,T],'j',R,[0,0,0,0,0,0,0,0,0,0,w,k]);
    end;
    wj = cell2mat(wj);
    wx = wj(:,1:N);
    wii = -0.5*wj(:,N+1:end)'*ainv(wx+wx')*wj(:,N+1:end);

    OFFLINE_TIME_DIST = toc
    RESIDUAL_1 = norm(WJ{1}-wx)
    RESIDUAL_2 = norm(WJ{2}-wii)

%% REDUCED ORDER MODEL EVALUATION
    l1 = zeros(1,N-1);
    l2 = zeros(1,N-1);
    l8 = zeros(1,N-1);

    for n=1:N-1
        uu = UU(:,1:n);
        y = ODE(LIN,OUT,[h,T],X,U,uu*uu'*P);
        l1(n) = norm(Y(:)-y(:),1)/n1;
        l2(n) = norm(Y(:)-y(:),2)/n2;
        l8(n) = norm(Y(:)-y(:),Inf)/n8;
    end;

%% PLOT REDUCED ORDER VS RELATIVE ERRORS
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

