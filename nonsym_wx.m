function nonsym_wx(o)
% nonsym_wx (nonsymmetric cross gramian linear state reduction)
% by Christian Himpe, 2013-2015 ( http://gramian.de )
% released under BSD 2-Clause License ( opensource.org/licenses/BSD-2-Clause )
%*
    if(exist('emgr')~=2)
        disp('emgr framework is required. Download at http://gramian.de/emgr.m');
        return;
    else
        global ODE;
        fprintf('emgr (version: %g)\n',emgr('version'));
    end

    %% Setup
    J = 8;
    N = 64;
    O = J;
    T = [0.0,0.01,1.0];
    L = (T(3)-T(1))/T(2);
    U = [ones(J,1),zeros(J,L-1)];
    X = zeros(N,1);

    rand('seed',1009);
    A = rand(N,N); A(1:N+1:end) = -0.55*N;
    B = rand(N,J);
    C = rand(O,N);

    LIN = @(x,u,p) A*x + B*u;
    OUT = @(x,u,p) C*x;

    norm1 = @(y) T(2)*sum(abs(y(:)));
    norm2 = @(y) sqrt(T(2)*dot(y(:),y(:)));
    norm8 = @(y) max(y(:));

    %% Main
    Y = ODE(LIN,OUT,T,X,U,0); % Full Order

    tic;
    WX = emgr(LIN,OUT,[J,N,O],T,'x',[],[0,0,0,0,0,0,1,0,0,0]);
    [UU,D,VV] = svd(WX);
    OFFLINE = toc

    for I=1:N-1
        uu = UU(:,1:I);
        vv = uu';
        a = vv*A*uu;
        b = vv*B;
        c = C*uu;
        x = vv*X;
        lin = @(x,u,p) a*x + b*u;
        out = @(x,u,p) c*x;
        y = ODE(lin,out,T,x,U,0); % Reduced Order
        l1(I) = norm1(Y-y)/norm1(Y);
        l2(I) = norm2(Y-y)/norm1(Y);
        l8(I) = norm8(Y-y)/norm1(Y);
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
