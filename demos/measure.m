function measure(o)
%%% summary: measure (input, state, output nonlinearity measure)
%%% project: emgr - Empirical Gramian Framework ( http://gramian.de )
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
    M = 1;		% number of inputs
    N = 8;		% number of states
    Q = M;		% number of outputs
    h = 0.01;		% time step size
    T = 1.0;		% time horizon

    rand('seed',1009);
    A = rand(N,N);		% \
    A(1:N+1:end) = -0.55*N;	%  system matrix 
    A = 0.5*(A+A');		% /
    B = rand(N,M);		% input matrix
    C = B';			% output matrix

    LIN = @(x,u,p,t) A*x + B*u;			% linear vector field
    OUT = @(x,u,p,t) C*x;			% linear output functional

    NIN = @(x,u,p,t) A*x + B*asinh(p*u);	% nonlinear input vector field
    NST = @(x,u,p,t) A*asinh(p*x) + B*u;	% nonlinear state vector field
    NOU = @(x,u,p,t) C*asinh(p*x);		% nonlinear output functional

%% GRAMIAN COMPUTATION
    tic;
    WL = emgr(LIN,OUT,[M,N,Q],[h,T],'x'); % linear reference gramian

    K = 20;
    y = zeros(3,K);
    P = 0;
    p = 2.0/K;

    for k=1:K
        Wi = emgr(NIN,OUT,[M,N,Q],[h,T],'c',P); % input nonlinearity gramian
        Ws = emgr(NST,OUT,[M,N,Q],[h,T],'x',P); % state nonlinearity gramian
        Wo = emgr(LIN,NOU,[M,N,Q],[h,T],'o',P); % output nonlinearity gramian

        Ni = norm(WL-Wi,'fro');
        Ns = norm(WL-Ws,'fro');
        No = norm(WL-Wo,'fro');

        y(:,k) = [Ni;Ns;No];
        P = P + p;
    end
    OFFLINE_TIME = toc

    y = y./trace(WL);

%% PLOT PARAMETER RANGE VS NONLINEARITY MEASURE
    if(nargin>0 && o==0), return; end; 
    figure('Name',mfilename,'NumberTitle','off');
    semilogy(linspace(0,2,K),y(1,:),'r','linewidth',2); hold on;
    semilogy(linspace(0,2,K),y(2,:),'g','linewidth',2);
    semilogy(linspace(0,2,K),y(3,:),'b','linewidth',2); hold off;
    pbaspect([2,1,1]);
    legend('Input ','State ','Output ','location','southeast');
    if(nargin>0 && o==1), print('-dsvg',[mfilename(),'.svg']); end;
end

