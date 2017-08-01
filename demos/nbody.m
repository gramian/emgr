function nbody(o)
%%% summary: nbody (5-body figure eight reduction)
%%% project: emgr - EMpirical Gramian Framework ( http://gramian.de )
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
    N = 5;		% number of particles
    h = 0.01;		% time step size
    T = 1.0;		% time horizon (offline)
    R = 4;		% target reduced number of particles
    p = ones(N,1);	% particle mass parameter

    X = [1.449;  0.0;    0.400; -0.345; -1.125;...
         0.448; -1.125; -0.448;  0.400;  0.345;...
         0.0;   -0.922; -1.335;  0.810; -0.919;...
        -0.349;  0.919; -0.349;  1.335;  0.810]; % 5-body eight-figure x0
    F = @(x,u,p,t) [x((2*N)+1:end);acc(x(1:2*N),u,p)];	% vector field
    G = @(x,u,p,t)  x(1:2*N);				% output functional

%% STRUCTURED REDUCED ORDER MODEL PROJECTION ASSEMBLY
    global ODE;
    ODE = @leapfrog;

    tic;
    WO = emgr(F,G,[0,4*N,2*N],[h,T],'o',p);
    WOP = WO(1:(2*N),1:(2*N));
    WOV = WO((2*N)+1:end,(2*N)+1:end);
    [PP,DD,QQ] = svd(WOP); PP = PP(:,1:2*R); QQ = PP'; %diag(DD)'
    [TT,DD,VV] = svd(WOV); TT = TT(:,1:2*R); VV = TT'; %diag(DD)'
    OFFLINE_TIME = toc

    T = 3.2;		% time horizon (online)
    U = @(t) 0;		% zero input function

    % Position-Based Reduction 
    f = @(x,u,p,t) [x((2*R)+1:end);QQ*acc(PP*x(1:2*R),u,p)];
    g = @(x,u,p,t) PP*x(1:2*R);
    x = [QQ*X(1:2*N);QQ*X((2*N)+1:end)];

    %{
    % Velocity-Based Reduction
    f = @(x,u,p) [x((2*R)+1:end);VV*acc(TT*x(1:2*R),u,p)];
    g = @(x,u,p) TT*x(1:2*R);
    x = [VV*X(1:2*N);VV*X((2*N)+1:end)];
    %}

%% REDUCED ORDER MODEL SIMULATION
    y = ODE(f,g,[h,T],x,U,p);
    ODE = [];

%% PLOT REDUCED ORDER MODEL TRAJECTORY
    if(nargin>0 && o==0), return; end; 
    figure('Name',mfilename,'NumberTitle','off');
    hold on;
    cmap = hsv(N+1);
    d = 1;
    for c=1:N
        plot(y(d,:),y(d+1,:),'--','Color',cmap(c,:),'linewidth',2);
        plot(y(d,end),y(d+1,end),'*','Color',cmap(c,:),'linewidth',2);
        d = d + 2;
    end
    hold off;
    ylim([-0.6,0.6]);
    pbaspect([2,1,1]);
    set(gcf,'InvertHardcopy','off');
    if(nargin>0 && o==1), print('-dsvg',[mfilename(),'.svg']); end;
end

function y = acc(x,u,p,t)
%%% summary: acc (acceleration vector-field)
%$
    N = numel(x)/2;
    A = reshape(x,[2,N]);
    y = zeros(2,N);

    for n=1:N
        B = bsxfun(@minus,A,A(:,n));
        Z = p'./(sqrt(1e-6+(sum(B.^2))).^3);
        B = bsxfun(@times,B,Z);
        y(:,n) = sum(B,2);
    end

    y = y(:);
end

function y = leapfrog(f,g,t,x0,u,p)
%%% summary: leapfrog (2nd order symplectic integrator)
%$
    h = t(1);
    L = floor(t(2)/h) + 1;
    N = numel(x0);
    n = N/2;

    y(:,1) = g(x0,u(0),p,0);
    y(end,L) = 0; % preallocate trajectory

    p1 = x0(1:n);
    q1 = x0(n+1:end);

    for l=2:L

        tl = (l-1)*h;
        ut = u(tl);
        a1 = f([p1;q1],ut,p,tl);
        p1 = p1 + h*q1 + 0.5*h*h*a1(n+1:end);
        a2 = f([p1;q1],ut,p,tl);
        q1 = q1 + 0.5*h*(a1(n+1:end) + a2(n+1:end));

        y(:,l) = g([p1;q1],ut,p,tl);
    end;
end
