function measure(o)
%%% summary: measure (input, state, output nonlinearity measure)
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
    N = 8;
    Q = M;
    T = [0.01,1.0];
    L = floor(T(2)/T(1)) + 1;
    X = zeros(N,1);

    rand('seed',1009);
    A = rand(N,N);
    A(1:N+1:end) = -0.55*N;
    A = 0.5*(A+A');
    B = rand(N,M);
    C = B';

    LIN = @(x,u,p,t) A*x + B*u;
    OUT = @(x,u,p,t) C*x;

    NIN = @(x,u,p,t) A*x + B*asinh(p*u);
    NST = @(x,u,p,t) A*asinh(p*x) + B*u;
    NOU = @(x,u,p,t) C*asinh(p*x);

%% OFFLINE
    tic;
    WL = emgr(LIN,OUT,[M,N,Q],T,'x'); % Linear Reference

    K = 20;
    y = zeros(3,K);
    P = 0;
    p = 2.0/K;

    for k=1:K
        Wi = emgr(NIN,OUT,[M,N,Q],T,'c',P); % Input Nonlinearity
        Ws = emgr(NST,OUT,[M,N,Q],T,'x',P); % State Nonlinearity
        Wo = emgr(LIN,NOU,[M,N,Q],T,'o',P); % Output Nonlinearity

        Ni = norm(WL-Wi,'fro');
        Ns = norm(WL-Ws,'fro');
        No = norm(WL-Wo,'fro');

        y(:,k) = [Ni;Ns;No];
        P = P + p;
    end
    OFFLINE = toc

    y = y./trace(WL);

%% OUTPUT
    if(nargin>0 && o==0), return; end; 
    figure('Name',mfilename,'NumberTitle','off');
    semilogy(linspace(0,2,K),y(1,:),'r','linewidth',2); hold on;
    semilogy(linspace(0,2,K),y(2,:),'g','linewidth',2);
    semilogy(linspace(0,2,K),y(3,:),'b','linewidth',2); hold off;
    pbaspect([2,1,1]);
    legend('Input ','State ','Output ','location','southeast');
    if(nargin>0 && o==1), print('-dsvg',[mfilename(),'.svg']); end;
end

