function benchmark_non(o)
% benchmark (nonlinear rc ladder)
% by Christian Himpe, 2013-2015 ( http://gramian.de )
% released under BSD 2-Clause License ( opensource.org/licenses/BSD-2-Clause )
%*
    if(exist('emgr')~=2)
        error('emgr not found! Get emgr at: http://gramian.de');
    else
        global ODE;
        fprintf('emgr (version: %g)\n',emgr('version'));
    end

    %% Setup
    J = 1;
    N = 64;
    O = 1;
    T = [0.0,0.01,1.0];
    L = (T(3)-T(1))/T(2);
    U = ones(1,L);
    X = zeros(N,1);

    g = @(x) exp(x)+x-1.0;

    A1 = spdiags(ones(N-1,1),-1,N,N)-speye(N);
    A2 = spdiags([ones(N-1,1);0],0,N,N)-spdiags(ones(N,1),1,N,N);

    NON = @(x,u,p) g(A1*x)-g(A2*x) + [u;sparse(N-1,1)];
    OUT = @(x,u,p) x(1);

    norm1 = @(y) T(2)*sum(abs(y(:)));
    norm2 = @(y) sqrt(T(2)*dot(y(:),y(:)));
    norm8 = @(y) max(y(:));

    %% Main
    Y = ODE(NON,OUT,T,X,U,0); % Full Order

    tic;
    WX = emgr(NON,OUT,[J,N,O],T,'x');
    [UU,D,VV] = svd(WX);
    OFFLINE = toc

    for I=1:N-1
        uu = UU(:,1:I);
        vv = uu';
        x = vv*X;
        non = @(x,u,p) vv*NON(uu*x,u,p);
        out = @(x,u,p) OUT(uu*x,u,p);
        y = ODE(non,out,T,x,U,0); % Reduced Order
        l1(I) = norm1(Y-y)/norm1(Y);
        l2(I) = norm2(Y-y)/norm2(Y);
        l8(I) = norm8(Y-y)/norm8(Y);
    end;

    %% Output
    if(nargin==0), return; end
    figure();
    semilogy(1:N-1,l1,'r','linewidth',2); hold on;
    semilogy(1:N-1,l2,'g','linewidth',2);
    semilogy(1:N-1,l8,'b','linewidth',2); hold off;
    xlim([1,N-1]);
    ylim([10^floor(log10(min([l1(:);l2(:);l8(:)]))-1),1]);
    pbaspect([2,1,1]);
    legend('L1 Error ','L2 Error ','L8 Error ','location','northeast');
    if(o==1), print('-dsvg',[mfilename(),'.svg']); end;
end
