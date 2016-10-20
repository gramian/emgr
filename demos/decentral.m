function decentral(o)
%%% summary: decentral (decentralized control)
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
    M = 4;
    N = 16;
    Q = M;
    T = [0.01,1.0];
    L = floor(T(2)/T(1)) + 1;
    U = @(t) ones(M,1)*(t<=T(1))/T(1);
    X = ones(N,1);

    rand('seed',1009);
    A = rand(N,N);
    A(1:N+1:end) = -0.55*N;
    B = rand(N,M);
    C = rand(Q,N);

    LIN = @(x,u,p,t) A*x + B*u;
    OUT = @(x,u,p,t) C*x;

%% OFFLINE

    global DOT;
    DOT = @(x,y) sum(sum(x.*y')); % Trace pseudo-kernel

    tic;
    PM = zeros(M,Q); % Participation Matrix

    for V=1:M
        for W=1:Q
            lin = @(x,u,p,t) A*x + B(:,V)*u;
            out = @(x,u,p,t) C(W,:)*x;
            io = emgr(lin,out,[1,N,1],T,'x');
            PM(V,W) = io(1,1);
        end
    end

    PM = PM./sum(sum(PM));
    OFFLINE = toc
    DOT = [];

%% OUTPUT
    if(nargin>0 && o==0), return; end; 
    figure('Name',mfilename,'NumberTitle','off');
    imagesc(PM);
    if(exist('viridis')==0), colormap(hot); end;
    colorbar;
    set(gca,'XTick',1:1:M,'YTick',1:1:Q);
    if(nargin>0 && o==1), print('-dpng',[mfilename(),'.png']); end;
end

