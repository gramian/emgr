function hierarchy(o)
% hierarchy (hierarchical network reduction)
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
    D = 3;	%Tree depth
    M = 4;	%Children per node

    J = 1;
    O = M^D;
    N = (M^(D+1)-1)/(M-1);
    T = [1,100];
    L = floor(T(2)/T(1)) + 1;
    U = exp(-0.0005*(1:L).^2);
    X = zeros(N,1);

    rand('seed',1009);
    A = trasm(D,M);
    B = sparse(N,1); B(1,1) = D;
    C = [sparse(O,N-O),speye(O)];

    LIN = @(x,u,p) A*x + B*u;
    OUT = @(x,u,p) C*x;

    ADJ = @(x,u,p) A'*x + C'*u;
    AOU = @(x,u,p) B'*x;

%% FULL ORDER
    ODE = @rk1;
    Y = rk1(LIN,OUT,T,X,U,0);
    n1 = norm(Y(:),1);
    n2 = norm(Y(:),2);
    n8 = norm(Y(:),Inf);

%% OFFLINE
    tic;
    WC = emgr(LIN,OUT,[J,N,O],T,'c');
    WO = emgr(ADJ,AOU,[O,N,J],T,'c');
    [VV,D,UU] = balance(WC,WO);
    OFFLINE = toc

%% EVALUATION
    for I=1:N-1
        uu = UU(:,1:I);
        vv = VV(1:I,:);
        a = vv*A*uu;
        b = vv*B;
        c = C*uu;
        x = vv*X;
        lin = @(x,u,p) a*x + b*u;
        out = @(x,u,p) c*x;
        y = ODE(lin,out,T,x,U,0);
        l1(I) = norm(Y(:)-y(:),1)/n1;
        l2(I) = norm(Y(:)-y(:),2)/n2;
        l8(I) = norm(Y(:)-y(:),Inf)/n8;
    end;

%% OUTPUT
    if(nargin>0 && o==0), return; end; 
    figure('Name',mfilename,'NumberTitle','off');
    semilogy(1:N-1,[l1;l2;l8],{'r','g','b'},'linewidth',2);
    xlim([1,N-1]);
    ylim([10^floor(log10(min([l1(:);l2(:);l8(:)]))-1),1]);
    pbaspect([2,1,1]);
    legend('L1 Error ','L2 Error ','L8 Error ','location','northeast');
    if(nargin>0 && o==1), print('-dsvg',[mfilename(),'.svg']); end;
end

%% ======== Tree Assembler ========
function A = trasm(d,c)

    a = (c^(d+1)-1)/(c-1);
    A = -speye(a);

    for I=0:((a-1)/c)-1
        b = 1+c*I;
        %A(1+I,b+1:b+c) = rand(1,c)+1;
        A(b+1:b+c,1+I) = rand(c,1)+1;
    end
end

%% ======== Balancer ========
function [X,Y,Z] = balance(WC,WO)

    [L,D,l] = svd(WC); LC = L*diag(sqrt(diag(D)));
    [L,D,l] = svd(WO); LO = L*diag(sqrt(diag(D)));
    [U,Y,V] = svd(LO'*LC);
    X = ( LO*U*diag(1.0./sqrt(diag(Y))) )';
    Z =   LC*V*diag(1.0./sqrt(diag(Y)));
end

%% ======== Integrator ========
function x = rk1(f,g,t,z,u,p)

    if(isnumeric(g) && g==1), g = @(x,u,p) x; end;

    h = t(1);
    L = floor(t(2)/h) + 1;

    x(:,1) = g(z,u(:,end),p);
    x(end,L) = 0;

    for l=2:L
        z = z + h*f(z,u(:,l-1),p);
        x(:,l) = g(z,u(:,l-1),p);
    end;
end

