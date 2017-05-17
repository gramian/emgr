function combined_wj(o)
%%% summary: combined_wj (joint gramian nonlinear combined reduction)
%%% project: emgr - Empirical Gramian Framework ( http://gramian.de )
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
    N = 64;				% number of states
    Q = M;				% number of outputs
    h = 0.01;				% time step size
    T = 1.0;				% time horizon
    X = zeros(N,1);			% initial state
    U = @(t) ones(M,1)*(t<=h)/h;	% impulse input function
    P = 0.9*rand(N,10) + 0.1;		% test parameters
    R = [0.1*ones(N,1),ones(N,1)];	% parameter range

    rand('seed',1009);
    A = full(sprand(N,N,1./N));		% \
    A(1:N+1:end) = -2.0;		%  system matrix
    A = A * 10.0;			% /
    B = 10.0*full(sprand(N,M,2.0/N));	% input matrix
    C = rand(Q,N);			% output matrix

    NON = @(x,u,p,t) A*tanh(p.*x) + B*u;	% vector field
    OUT = @(x,u,p,t) C*x;			% output functional

%% REDUCED ORDER MODEL PROJECTION ASSEMBLY
    tic;
    WJ = emgr(NON,OUT,[M,N,Q],[h,T],'j',R,[0,0,0,0,0,0,0,1,1,0,0,0]);
    [UU,D,VV] = svd(WJ{1});
    [PP,D,QQ] = svd(WJ{2});
    OFFLINE = toc

%% REDUCED ORDER MODEL EVALUATION
    Q = size(P,2);
    s = 4;
    l2 = zeros(N/s,N/s,Q);

    for q=1:Q
        Y = ODE(NON,OUT,[h,T],X,U,P(:,q));
        n2 = norm(Y(:),2);

        a = 1;
        for n=1:s:N
            uu = UU(:,1:n);
            vv = uu';
            b = 1;
            for k=1:s:N
                pp = PP(:,1:k);
                qq = pp';
                p = pp*qq*P(:,q);
                x = vv*X;
                non = @(x,u,p,t) vv*NON(uu*x,u,p);
                out = @(x,u,p,t) OUT(uu*x,u,p);
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
    xlim([1,N/s]);
    ylim([1,N/s]);
    set(gca,'ZScale','log');
    zlim([1e-16,1]);
    set(gca,'XTick',1:2:N/s,'XTickLabel',1:N/(2*s):N);
    set(gca,'YTick',1:2:N/s,'YTickLabel',1:N/(2*s):N);
    ylabel('State Dimension')
    xlabel('Parameter Dimension');
    set(h,'CData',log10(get(h,'CData')))
    view(135,30);
    if(nargin>0 && o==1), print('-dsvg',[mfilename(),'.svg']); end;
end

