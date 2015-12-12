function benchmark_lin(o)
% benchmark (iss model reduction benchmark)
% by Christian Himpe, 2013-2015 ( http://gramian.de )
% released under BSD 2-Clause License ( opensource.org/licenses/BSD-2-Clause )
%*
    if(exist('emgr')~=2)
        error('emgr not found! Get emgr at: http://gramian.de');
    else
        global ODE; ODE = [];
        fprintf('emgr (version: %g)\n',emgr('version'));
    end

    D = 'iss';
    if(exist(['/tmp/',D,'.mat'],'file')==0)
        urlwrite(['http://slicot.org/objects/software/shared/bench-data/',D,'.zip'],['/tmp/',D,'.zip']);
        unzip(['/tmp/',D,'.zip'],'/tmp');
    end
    load(['/tmp/',D,'.mat']);

%% SETUP
    N = size(A,2);
    n = N/2;
    J = size(B,2);
    O = size(C,1);
    T = [0.01,1.0];
    L = floor(T(2)/T(1)) + 1;
    U = [ones(J,1),zeros(J,L-1)];
    X = zeros(N,1);

    LIN = @(x,u,p) A*x + B*u;
    OUT = @(x,u,p) C*x;

%% FULL ORDER
    Y = ODE(LIN,OUT,T,X,U,0);
    n1 = norm(Y(:),1);
    n2 = norm(Y(:),2);
    n8 = norm(Y(:),Inf);

%% OFFLINE
    tic;
    WX = emgr(LIN,OUT,[J,N,O],T,'x');
    %WXP = WX(1:n,1:n);
    %[UU D VV] = svd(WXP);
    WXV = WX(n+1:N,n+1:N);
    [UU,D,VV] = svd(WXV);
    OFFLINE = toc

%% EVALUATION
    for I=1:n-1
        uu = UU(:,1:I);
        vv = uu';
        a = [zeros(I,I),eye(I);vv*A(n+1:N,1:n)*uu,vv*A(n+1:N,n+1:N)*uu];
        b = [zeros(I,J);vv*B(n+1:N,:)];
        c = [zeros(O,I),C(:,n+1:N)*uu];
        x = zeros(2*I,1);
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
    semilogy(2:2:N-2,l1,'r','linewidth',2); hold on;
    semilogy(2:2:N-2,l2,'g','linewidth',2);
    semilogy(2:2:N-2,l8,'b','linewidth',2); hold off;
    xlim([2,N-2]);
    ylim([10^floor(log10(min([l1(:);l2(:);l8(:)]))-1),1]);
    pbaspect([2,1,1]);
    legend('L1 Error ','L2 Error ','L8 Error ','location','northeast');
    if(nargin>0 && o==1), print('-dsvg',[mfilename(),'.svg']); end;
end
