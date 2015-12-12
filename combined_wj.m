function combined_wj(o)
% combined_wj (joint gramian nonlinear combined reduction)
% by Christian Himpe, 2013-2015 ( http://gramian.de )
% released under BSD 2-Clause License ( opensource.org/licenses/BSD-2-Clause )
%*
    if(exist('emgr')~=2)
        error('emgr not found! Get emgr at: http://gramian.de');
    else
        global ODE; ODE = [];
        fprintf('emgr (version: %g)\n',emgr('version'));
    end

%% SETUP
    J = 1;
    N = 64;
    O = J;
    T = [0.01,1.0];
    L = floor(T(2)/T(1)) + 1;
    U = [ones(J,1)./T(2),zeros(J,L)];
    X = zeros(N,1);

    rand('seed',1009);

    A = rand(N,N);
    A(1:N+1:end) = -0.9*N; 
    B = 30.0*full(sprand(N,J,0.1));
    C = rand(O,N);

    P = 0.5*rand(N,1) + 0.5;
    Q = [0.5*ones(N,1),ones(N,1)];

    NON = @(x,u,p) A*tanh(p.*x) + B*u;
    OUT = @(x,u,p) C*x;

%% FULL ORDER
    Y = ODE(NON,OUT,T,X,U,P);
    n2 = norm(Y(:),2);

%% OFFLINE
    tic;
    WJ = emgr(NON,OUT,[J,N,O],T,'j',Q,[0,0,0,0,0,0,0,0,1,1,1,0]);
    [UU,D,VV] = svd(WJ{1});
    [PP,D,QQ] = svd(WJ{2});
    OFFLINE = toc

%% EVALUATION
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
            y = ODE(non,out,T,x,U,p);
            l2(a,b) = min(1.0,norm(Y(:)-y(:),2)/n2);
            b = b + 1;
        end;
        a = a + 1;
    end;

%% OUTPUT
    if(nargin>0 && o==0), return; end; 
    figure('Name',mfilename,'NumberTitle','off');
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
    if(nargin>0 && o==1), print('-dsvg',[mfilename(),'.svg']); end;
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
