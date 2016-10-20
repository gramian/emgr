function nbody(o)
%%% summary: nbody (5-body figure eight reduction)
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
    N = 5;
    T = [0.01,1.0];
    R = 4;
    p = ones(N,1);

    X = [1.449;  0.0;    0.400; -0.345; -1.125;...
         0.448; -1.125; -0.448;  0.400;  0.345;...
         0.0;   -0.922; -1.335;  0.810; -0.919;...
        -0.349;  0.919; -0.349;  1.335;  0.810]; %5-body eight-figure
    F = @(x,u,p,t) [x((2*N)+1:end);acc(x(1:2*N),u,p)];
    G = @(x,u,p,t)  x(1:2*N);

%% OFFLINE
    tic;
    WO = emgr(F,G,[0,4*N,2*N],T,'o',p);
    WOP = WO(1:(2*N),1:(2*N));
    WOV = WO((2*N)+1:end,(2*N)+1:end);
    [PP,DD,QQ] = svd(WOP); PP = PP(:,1:2*R); QQ = PP'; %diag(DD)'
    [TT,DD,VV] = svd(WOV); TT = TT(:,1:2*R); VV = TT'; %diag(DD)'
    OFFLINE = toc

    T = [0.01,3.0];
    L = floor(T(2)/T(1)) + 1;
    U = @(t) 0;

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

%% ONLINE
    y = ODE(f,g,T,x,U,p);

%% OUTPUT
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

function y = acc(x,u,p)
%%% summary: acc (acceleration vector-field)
%$
    N = numel(x)/2;
    A = reshape(x,[2,N]);
    B = zeros(size(A));
    Z = zeros(1,N);
    y = zeros(2,N);

    for n=1:N
        B = bsxfun(@minus,A,A(:,n));
        Z = p'./(sqrt(1e-6+(sum(B.^2))).^3);
        B = bsxfun(@times,B,Z);
        y(:,n) = sum(B,2);
    end

    y = y(:);
end

