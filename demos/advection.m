function advection(o)
%%% summary: advection (finite difference discretized transport equation)
%%% project: emgr - Empirical Gramian Framework ( gramian.de )
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
    M = 0;				% number of inputs
    N = 256;				% number of states
    Q = N;				% number of outputs
    R = 10;				% target reduced order
    h = 0.01;				% time step size
    T = 0.1;				% time horizon
    X = exp(-linspace(-2,8,N).^2)';	% initial state
    U = @(t) 0;				% zero input function
    P = 0.55;				% parameter

    A = spdiags(N*[ones(N,1),-ones(N,1)],[-1,0],N,N);	% system matrix

    LIN = @(x,u,p,t) p*A*x;		% vector field
    ADJ = @(x,u,p,t) p*A'*x + u;	% adjoint vector field
    OUT = @(x,u,p,t) x;			% output functional

%% REDUCED ORDER MODEL PROJECTION ASSEMBLY
    tic;
    WO = emgr(ADJ,1,[N,N,Q],[h,T],'c',P);
    [UU,D,VV] = svd(WO);
    UU = UU(:,1:R);
    VV = UU';
    a = VV*A*UU;
    c = UU;
    x = VV*X;
    lin = @(x,u,p,t) p*a*x;
    out = @(x,u,p,t) c*x;
    OFFLINE = toc

%% REDUCED ORDER MODEL SIMULATION
    y = ODE(lin,out,[0.01,1.0],x,U,P);

%% PLOT REDUCED ORDER MODEL PHASE SPACE
    if(nargin>0 && o==0), return; end; 
    figure('Name',mfilename,'NumberTitle','off');
    imagesc(sparse(y)); caxis([0,max(y(:))]);
    set(gca,'YTick',0,'xtick',[]); ylabel('X'); xlabel('t');
    pbaspect([2,1,1]);
    if(nargin>0 && o==1), print('-dpng',[mfilename(),'.png']); end;
end

