function test_cws(o)
%%% summary: test_cws (sensitivity gramian combined reduction)
%%% project: emgr - Empirical Gramian Framework ( http://gramian.de )
%%% authors: Christian Himpe ( 0000-0003-2194-6754 )
%%% license: 2-Clause BSD (2016--2017)
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
    r = 0.5*ones(N,1);			% nominal parameter

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

%% REDUCED ORDER MODEL PROJECTION ASSEMBLY
    tic;
    WS = emgr(LIN,OUT,[M,N,Q],[h,T],'s',R);    
    WO = emgr(LIN,OUT,[M,N,Q],[h,T],'o',r);
    [L1,D1,R1] = svd(WS{1}); LC = L1*diag(sqrt(diag(D1)));
    [L2,D2,R2] = svd(WO); LO = L2*diag(sqrt(diag(D2)));
    [Lb,Db,Rb] = svd(LO'*LC);
    UU = ( LO*Lb*diag(1.0./sqrt(diag(Db))) )';
    VV =   LC*Rb*diag(1.0./sqrt(diag(Db)));
    [D,V] = sort(WS{2},'descend');
    PP = eye(N);
    PP = PP(:,V);
    OFFLINE = toc

%% REDUCED ORDER MODEL EVALUATION
    l1 = zeros(1,N-1);
    l2 = zeros(1,N-1);
    l8 = zeros(1,N-1);

    for n=1:N-1
        uu = VV(:,1:n);
        vv = UU(1:n,:);
        pp = PP(:,1:n);
        lin = @(x,u,p,t) vv*LIN(uu*x,u,p,t);
        out = @(x,u,p,t) OUT(uu*x,u,p,t);
        y = ODE(lin,out,[h,T],vv*X,U,pp*pp'*P);
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
    legend('L1 Error ','L2 Error ','L8 Error ','location','southeast');
    if(nargin>0 && o==1), print('-dsvg',[mfilename(),'.svg']); end;
end

