function benchmark_non(o)
%%% summary: benchmark_non (nonlinear rc ladder benchmark)
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

%% SYSTEM SETUP
    M = 1;				% number of inputs
    N = 64;				% number of states
    Q = M;				% number of outputs
    h = 0.01;				% time step size
    T = 1.0;				% time horizon
    X = zeros(N,1);			% initial state
    U = @(t) ones(M,1)*(t<=0.5*T);	% rectangle input function

    g = @(x) exp(x) + x - 1.0;		% diode nonlinearity

    A0 = sparse(N,N); A0(1,1) = 1;
    A1 = spdiags(ones(N-1,1),-1,N,N)-speye(N); A1(1,1) = 0;
    A2 = spdiags([ones(N-1,1);0],0,N,N)-spdiags(ones(N,1),1,N,N);
    B = sparse(N,1); B(1,1) = 1;

    NON = @(x,u,p,t) -g(A0*x) + g(A1*x) - g(A2*x) + B*u; % vector field
    OUT = @(x,u,p,t) x(1);				 % output functional

%% FULL ORDER MODEL REFERENCE SOLUTION
    Y = ODE(NON,OUT,[h,T],X,U,0);
    n1 = norm(Y(:),1);
    n2 = norm(Y(:),2);
    n8 = norm(Y(:),Inf);

%% REDUCED ORDER MODEL PROJECTION ASSEMBLY
    tic;
    WX = emgr(NON,OUT,[M,N,Q],[h,T],'x',0,0,U);
    [UU,D,VV] = svd(WX);
    OFFLINE_TIME = toc

%% REDUCED ORDER MODEL EVALUATION
    l1 = zeros(1,N-1);
    l2 = zeros(1,N-1);
    l8 = zeros(1,N-1);

    for n=1:N-1
        uu = UU(:,1:n);
        vv = uu';
        x = vv*X;
        non = @(x,u,p,t) vv*NON(uu*x,u,p);
        out = @(x,u,p,t) OUT(uu*x,u,p);
        y = ODE(non,out,[h,T],x,U,0);
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

