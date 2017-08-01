function advection(o)
%%% summary: advection (finite difference discretized transport equation)
%%% project: emgr - EMpirical GRamian Framework ( http://gramian.de )
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

%% SYSTEM SETUP
    M = 1;					% number of inputs
    N = 256;					% number of states
    Q = 1;					% number of outputs
    L = N;					% number of time steps
    T = 1.0;					% time horizon
    h = T./L;					% time step width
    p = 1.3;					% velocity
    A = N * spdiags([ones(N,1),-ones(N,1)],[-1,0],N,N); % discretization
    B = [1.0;sparse(N-1,1)];			% input matrix
    C = [sparse(1,N-1),1.0];			% output matrix
    X = zeros(N,1);				% initial state
    U = @(t) exp(((t-0.1).^2)./(-0.001));	% input function

    LIN = @(x,u,p,t) p*A*x + B*u;		% vector field
    ADJ = @(x,u,p,t) p*A'*x + C'*u;		% adjoint vector field
    OUT = @(x,u,p,t) C*x;			% output functional

%% FULL ORDER MODEL REFERENCE SOLUTION
    Y = ODE(LIN,OUT,[h,T],X,U,p);
    n1 = norm(Y(:),1);
    n2 = norm(Y(:),2);
    n8 = norm(Y(:),Inf);

%% REDUCED ORDER MODEL PROJECTION ASSEMBLY
    tic;
    WX = emgr(LIN,ADJ,[M,N,Q],[h,T],'y',p,[4,0,0,0,0,0,0,0,0,0,0,0]);
    [UU,D,VV] = svd(WX);
    OFFLINE_TIME = toc

%% REDUCED ORDER MODEL EVALUATION
    l1 = zeros(1,N-1);
    l2 = zeros(1,N-1);
    l8 = zeros(1,N-1);

    for n=1:4:N-1
        uu = UU(:,1:n);
        a = uu'*A*uu;
        b = uu'*B;
        c = C*uu;
        lin = @(x,u,p,t) p*a*x + b*u;
        out = @(x,u,p,t) c*x;
        y = ODE(lin,out,[h,T],uu'*X,U,p);
        l1(n) = norm(Y(:)-y(:),1)/n1;
        l2(n) = norm(Y(:)-y(:),2)/n2;
        l8(n) = norm(Y(:)-y(:),Inf)/n8;
    end;

%% PLOT REDUCDED ORDER VS RELATIVE ERRORS
    if(nargin>0 && o==0), return; end; 
    figure('Name',mfilename,'NumberTitle','off');
    semilogy(1:4:N-1,l1(1:4:N-1),'r','linewidth',2); hold on;
    semilogy(1:4:N-1,l2(1:4:N-1),'g','linewidth',2);
    semilogy(1:4:N-1,l8(1:4:N-1),'b','linewidth',2); hold off;
    xlim([1,N-1]);
    ylim([1e-16,1]);
    pbaspect([2,1,1]);
    legend('L1 Error ','L2 Error ','L8 Error ','location','northeast');
    if(nargin>0 && o==1), print('-dsvg',[mfilename(),'.svg']); end;
end

