function measure(o)
% measure (nonlinearity measure)
% by Christian Himpe, 2013-2016 ( http://gramian.de )
% released under BSD 2-Clause License ( opensource.org/licenses/BSD-2-Clause )
%*
    if(exist('emgr')~=2)
        error('emgr not found! Get emgr at: http://gramian.de');
    else
        global ODE; ODE = [];
        fprintf('emgr (version: %g)\n',emgr('version'));
    end

%% SETUP
    J = 1;
    N = 8;
    O = J;
    T = [0.01,1.0];
    L = floor(T(2)/T(1)) + 1;
    X = zeros(N,1);

    rand('seed',1009);
    A = rand(N,N);
    A(1:N+1:end) = -0.55*N;
    A = 0.5*(A+A');
    B = rand(N,J);
    C = B';

    LIN = @(x,u,p) A*x + B*u;
    OUT = @(x,u,p) C*x;

    NIN = @(x,u,p) A*x + B*asinh(p*u);
    NST = @(x,u,p) A*asinh(p*x) + B*u;
    NOU = @(x,u,p) C*asinh(p*x);

%% OFFLINE
    tic;
    WL = emgr(LIN,OUT,[J,N,O],T,'x'); % Linear Reference

    K = 20;
    y = zeros(3,K);
    Q = 2.0/K;
    P = 0;

    for(I=1:K)
        Wi = emgr(NIN,OUT,[J,N,O],T,'c',P); % Input Nonlinearity
        Ws = emgr(NST,OUT,[J,N,O],T,'x',P); % State Nonlinearity
        Wo = emgr(LIN,NOU,[J,N,O],T,'o',P); % Output Nonlinearity

        Ni = sum(sum(abs(WL-Wi)));
        Ns = sum(sum(abs(WL-Ws)));
        No = sum(sum(abs(WL-Wo)));

        y(:,I) = [Ni;Ns;No];
        P = P + Q;
    end
    OFFLINE = toc

    y = y./trace(WL);

%% OUTPUT
    if(nargin>0 && o==0), return; end; 
    figure('Name',mfilename,'NumberTitle','off');
    semilogy(linspace(0,2,K),y,{'r','g','b'},'linewidth',2);
    pbaspect([2,1,1]);
    legend('Input ','State ','Output ','location','southeast');
    if(nargin>0 && o==1), print('-dsvg',[mfilename(),'.svg']); end;
end
