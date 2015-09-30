function benchmark_lin(o)
% benchmark (iss model reduction benchmark)
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

    D = 'iss';
    if(exist(['/tmp/',D,'.mat'],'file')==0)
        urlwrite(['http://slicot.org/objects/software/shared/bench-data/',D,'.zip'],['/tmp/',D,'.zip']);
        unzip(['/tmp/',D,'.zip'],'/tmp');
    end
    load(['/tmp/',D,'.mat']);

    %% Setup
    N = size(A,2);
    n = N/2;
    J = size(B,2);
    O = size(C,1);
    T = [0.0,0.01,1.0];
    L = (T(3)-T(1))/T(2);
    U = [ones(J,1),zeros(J,L-1)];
    X = zeros(N,1);

    LIN = @(x,u,p) A*x + B*u;
    OUT = @(x,u,p) C*x;

    norm1 = @(y) T(2)*sum(abs(y(:)));
    norm2 = @(y) sqrt(T(2)*dot(y(:),y(:)));
    norm8 = @(y) max(y(:));

    %% Main
    Y = ODE(LIN,OUT,T,X,U,0); % Full Order

    tic;
    WX = emgr(LIN,OUT,[J,N,O],T,'x');
    %WXP = WX(1:n,1:n);
    %[UU D VV] = svd(WXP);
    WXV = WX(n+1:N,n+1:N);
    [UU,D,VV] = svd(WXV);
    OFFLINE = toc

    for I=1:n-1
        uu = UU(:,1:I);
        vv = uu';
        a = [zeros(I,I),eye(I);vv*A(n+1:N,1:n)*uu,vv*A(n+1:N,n+1:N)*uu];
        b = [zeros(I,J);vv*B(n+1:N,:)];
        c = [zeros(O,I),C(:,n+1:N)*uu];
        x = zeros(2*I,1);
        lin = @(x,u,p) a*x + b*u;
        out = @(x,u,p) c*x;
        y = ODE(lin,out,T,x,U,0); % Reduced Order
        l1(I) = norm1(Y-y)/norm1(Y);
        l2(I) = norm2(Y-y)/norm2(Y);
        l8(I) = norm8(Y-y)/norm8(Y);
    end;

    %% Output
    if(nargin==0), return; end
    figure();
    semilogy(2:2:N-2,l1,'r','linewidth',2); hold on;
    semilogy(2:2:N-2,l2,'g','linewidth',2);
    semilogy(2:2:N-2,l8,'b','linewidth',2); hold off;
    xlim([2,N-2]);
    ylim([10^floor(log10(min([l1(:);l2(:);l8(:)]))-1),1]);
    pbaspect([2,1,1]);
    legend('L1 Error ','L2 Error ','L8 Error ','location','northeast');
    if(o==1), print('-dsvg',[mfilename(),'.svg']); end;
end
