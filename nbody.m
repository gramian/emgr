function nbody(o)
% nbody (n-body reduction)
% by Christian Himpe, 2013-2015 ( http://gramian.de )
% released under BSD 2-Clause License ( opensource.org/licenses/BSD-2-Clause )
%*
%*
    if(exist('emgr')~=2)
        error('emgr not found! Get emgr at: http://gramian.de');
    else
        global ODE; ODE = [];
        fprintf('emgr (version: %g)\n',emgr('version'));
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
    F = @(x,u,p) [x((2*N)+1:end);acc(x(1:2*N),u,p)];
    G = @(x,u,p)  x(1:2*N);

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
    U = sparse(1,L);

    % Position-Based Reduction 
    f = @(x,u,p) [x((2*R)+1:end);QQ*acc(PP*x(1:2*R),u,p)];
    g = @(x,u,p) PP*x(1:2*R);
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

%% ======== Acceleration ========
function y = acc(x,u,p)

    N = numel(x)/2;
    A = reshape(x,[2,N]);
    B = zeros(size(A));
    Z = zeros(1,N);
    y = zeros(2,N);

    for I=1:N
        B = bsxfun(@minus,A,A(:,I));
        Z = p'./(sqrt(1e-6+(sum(B.^2))).^3);
        B = bsxfun(@times,B,Z);
        y(:,I) = sum(B,2);
    end

    y = y(:);
end
