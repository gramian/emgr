function advection(o)
%%% summary: advection (finite difference discretized transport equation)
%%% project: emgr - Empirical Gramian Framework ( gramian.de )
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
    M = 0;
    N = 256;
    Q = N;
    R = 10;
    T = [0.01,0.1];
    L = floor(T(2)/T(1)) + 1;
    U = @(t) 0;
    X = exp(-linspace(-2,8,N).^2)';

    P = 0.55;
    A = spdiags(N*[ones(N,1),-ones(N,1)],[-1,0],N,N);

    LIN = @(x,u,p,t) p*A*x;
    ADJ = @(x,u,p,t) p*A'*x + u;
    OUT = @(x,u,p,t) x;

%% OFFLINE
    tic;
    WO = emgr(ADJ,1,[N,N,Q],T,'c',P);
    [UU,D,VV] = svd(WO); UU = UU(:,1:R); VV = UU';
    a = VV*A*UU;
    c = UU;
    x = VV*X;
    lin = @(x,u,p,t) p*a*x;
    out = @(x,u,p,t) c*x;
    OFFLINE = toc

%% ONLINE
    y = ODE(lin,out,[0.01,1.0],x,U,P);

%% OUTPUT
    if(nargin>0 && o==0), return; end; 
    figure('Name',mfilename,'NumberTitle','off');
    imagesc(sparse(y)); caxis([0,max(y(:))]);
    if(exist('viridis')==0), colormap(hot); end;
    set(gca,'YTick',0,'xtick',[]); ylabel('X'); xlabel('t');
    pbaspect([2,1,1]);
    if(nargin>0 && o==1), print('-dpng',[mfilename(),'.png']); end;
end

