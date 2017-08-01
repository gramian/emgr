function energy_wz(o)
%%% summary: energy_wz (cross gramian experimental reduction)
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
    M = 4;				% number of inputs
    N = M*M*M;				% number of states
    Q = 1;				% number of outputs
    h = 0.01;				% time step size
    T = 1.0;				% time horizon
    U = @(t) ones(M,1)*(t<=h)/h;	% impulse input function
    X = zeros(N,1);			% initial state

    rand('seed',1009);		% \
    A = rand(N,N);		%  system matrix
    A(1:N+1:end) = -0.55*N;	% /
    B = rand(N,M);		% input matrix
    C = rand(Q,N);		% output matrix

    LIN = @(x,u,p,t) A*x + B*u;	% vector field
    OUT = @(x,u,p,t) x'*x;	% output functional

%% FULL ORDER MODEL REFERENCE SOLUTION
    Y = ODE(LIN,OUT,[h,T],X,U,0);
    n1 = norm(Y(:),1);
    n2 = norm(Y(:),2);
    n8 = norm(Y(:),Inf);

% REDUCED ORDER MODEL PROJECTION ASSEMBLY
    tic;
    WZ = emgr(LIN,OUT,[M,N,Q],[h,T],'x',0,[0,0,0,0,0,0,1,0,0,0,0,0],1,0,1);
    [UU,D,VV] = svd(WZ);
    OFFLINE = toc

%% REDUCED ORDER MODEL EVALUATION
    l1 = zeros(1,N-1);
    l2 = zeros(1,N-1);
    l8 = zeros(1,N-1);

    for n=1:N-1
        uu = UU(:,1:n);
        vv = uu';
        a = vv*A*uu;
        b = vv*B;
        x = vv*X;
        lin = @(x,u,p,t) a*x + b*u;
        out = @(x,u,p,t) OUT(uu*x,u,p);
        y = ODE(lin,out,[h,T],x,U,0);
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

