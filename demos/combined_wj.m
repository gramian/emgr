function combined_wj(o)
%%% summary: combined_wj (joint gramian nonlinear combined reduction)
%%% project: emgr - EMpirical GRamian Framework ( http://gramian.de )
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
    M = 1;				% number of inputs
    N = 96;				% number of states
    Q = M;				% number of outputs
    h = 0.01;				% time step size
    T = 1.0;				% time horizon
    X = zeros(N,1);			% initial state
    U = @(t) ones(M,1)*(t<=h)/h;	% impulse input function
    P = 0.5*rand(N,10) + 0.5;		% test parameters
    R = [0.5*ones(N,1),ones(N,1)];	% parameter range

    rand('seed',1009);
    A = 0.1*N*toeplitz([-4,1,zeros(1,N-2)]);	% system matrix
    B = 1e-16 * ( 1e16.^rand(N,M) );		% input matrix
    C = ones(Q,N);				% output matrix

    NON = @(x,u,p,t) A*tanh(p.*x) + B*u;	% vector field
    OUT = @(x,u,p,t) C*x;			% output functional

%% REDUCED ORDER MODEL PROJECTION ASSEMBLY
    tic;
    WJ = emgr(NON,OUT,[M,N,Q],[h,T],'j',R,[3,0,0,0,0,0,0,1,1,0,0,0]);
    [UU,D,VV] = svd(WJ{1});
    [PP,D,QQ] = svd(WJ{2});
    OFFLINE_TIME = toc

%% REDUCED ORDER MODEL EVALUATION
    Q = size(P,2);
    s = 16;
    l2 = zeros(N/s+1,N/s+1,Q);

    for q=1:Q
        Y = ODE(NON,OUT,[h,T],X,U,P(:,q));
        n2 = norm(Y(:),2);

        a = 1;
        for n=[1,s:s:N]
            uu = UU(:,1:n);
            vv = uu';
            non = @(x,u,p,t) vv*NON(uu*x,u,p);
            out = @(x,u,p,t) OUT(uu*x,u,p);
            x = vv*X;

            b = 1;
            for k=[1,s:s:N]
                pp = PP(:,1:k);
                qq = pp';
                p = pp*qq*P(:,q);
                y = ODE(non,out,[h,T],x,U,p);
                l2(a,b,q) = norm(Y(:)-y(:),2) / n2;
                b = b + 1;
            end;
            a = a + 1;
        end;
    end;
    l2 = sqrt(sum(l2.^2,3));

%% PLOT REDUCED ORDER VS RELATIVE JOINT L2 ERROR
    if(nargin>0 && o==0), return; end; 
    figure('Name',mfilename,'NumberTitle','off');
    h = surf(min(1.0,l2));
    xlim([1,N/s+1]);
    ylim([1,N/s+1]);
    set(gca,'ZScale','log');
    zlim([1e-4,1]);
    set(gca,'XTick',[3:2:N/s],'XTickLabel',[2*s:2*s:N]);
    set(gca,'YTick',[3:2:N/s],'YTickLabel',[2*s:2*s:N]);
    ylabel('State Dimension')
    xlabel('Parameter Dimension');
    set(h,'CData',log10(get(h,'CData')));
    set(gca,'CLim',log10(get(gca,'ZLim')));
    view(135,30);
    if(nargin>0 && o==1), print('-dsvg',[mfilename(),'.svg']); end;
end

