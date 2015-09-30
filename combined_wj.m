function combined_wj(o)
% combined_wj (joint gramian nonlinear combined reduction)
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
    randn('seed',1009);
    A = rand(N,N); A(1:N+1:end) = -N; A = 0.5*(A+A');
    B = rand(N,J);
    C = B';
    P = 0.1*randn(N,1)+0.5;

    NON = @(x,u,p) A*tanh(p.*x) + B*u;
    OUT = @(x,u,p) C*x;

    norm2 = @(y) sqrt(T(2)*dot(y(:),y(:)));

    %% Main
    Y = ODE(NON,OUT,T,X,U,P);

    l2 = zeros(16,16);

    tic;
    WJ = emgr(NON,OUT,[J,N,O],T,'j',P,0,1,0,1);
    [UU,D,VV] = svd(WJ{1});
    [PP,D,QQ] = svd(WJ{2});
    OFFLINE = toc

    a = 1;
    for I=1:4:N
        uu = UU(:,1:I);
        vv = uu';
        b = 1;
        for K=1:4:N
            pp = PP(:,1:K);
            qq = pp';
            p = pp*qq*P;
            x = vv*X;
            non = @(x,u,p) vv*NON(uu*x,u,p);
            out = @(x,u,p) OUT(uu*x,u,p);
            y = ODE(non,out,T,x,U,p); % Reduced Order
            l2(a,b) = norm2(Y-y)/norm2(Y);
            b = b + 1;
        end;
        a = a + 1;
    end;

    %% Output
    if(nargin==0), return; end
    figure();
    h = surf(l2);
    xlim([1,16]);
    ylim([1,16]);
    set(gca,'ZScale','log');
    set(gca,'XTick',1:2:16,'XTickLabel',1:8:N);
    set(gca,'YTick',1:2:16,'YTickLabel',1:8:N);
    ylabel('State Dimension')
    xlabel('Parameter Dimension');
    colormap(antijet);
    set(h,'CData',log10(get(h,'CData')))
    view(135,30);
    if(o==1), print('-dsvg',[mfilename(),'.svg']); end;
end

%% ======== Colormap ========
function m = antijet(n)
% antijet colormap
% by Christian Himpe 2014-2015
% released under BSD 2-Clause License ( opensource.org/licenses/BSD-2-Clause )
%*
    if(nargin<1 || isempty(n)), n = 256; end;
    L = linspace(0,1,n);

    R = -0.5*sin( L*(1.37*pi)+0.13*pi )+0.5;
    G = -0.4*cos( L*(1.5*pi) )+0.4;
    B = 0.3*sin( L*(2.11*pi) )+0.3;

    m = [R;G;B]';
end
