function indices(o)
%%% summary: indices (cross-gramian-based system indices)
%%% project: emgr - EMpirical GRamian Framework ( http://gramian.de )
%%% authors: Christian Himpe ( 0000-0003-2194-6754 )
%%% license: 2-Clause BSD (2017)
%$
    if(exist('emgr')~=2)
        error('emgr not found! Get emgr at: http://gramian.de');
    else
        global ODE;
        ODE = [];
        fprintf('emgr (version: %1.1f)\n',emgr('version'));
    end

    rand('seed',1009);

%% SYSTEM SETUP
    M = 1;					% number of inputs
    N = 64;					% number of states
    Q = 1;					% number of outputs

    T = 1.0;					% time horizon
    h = 0.01;					% time step width
    A = 1.3*N * spdiags([ones(N,1),-ones(N,1)],[-1,0],N,N);			 % discretization
    B = rand(N,1);				% input matrix
    C = rand(1,N);				% output matrix
    X = zeros(N,1);				% initial state
    U = @(t) t<=h;				% input function

    LIN = @(x,u,p,t) A*x + B*u;			% vector field
    ADJ = @(x,u,p,t) A'*x + C'*u;		% adjoint vector field
    OUT = @(x,u,p,t) C*x;			% output functional

%figure; plot(0:h:T,ODE(LIN,OUT,[h,T],X,U,0)); return;

%% REDUCED ORDER MODEL PROJECTION ASSEMBLY
    tic;
    WX = emgr(LIN,ADJ,[M,N,Q],[h,T],'y');

    [EV,ev] = eig(WX);
    ev = diag(ev);
    [ev,pm] = sort(ev,'descend');
    EV = EV(:,pm);
    sv = svd(WX);
    b = EV\B;
    c = C*EV;
    OFFLINE_TIME = toc

%% REDUCED ORDER MODEL EVALUATION
    hn = zeros(1,N-1);
    l1 = zeros(1,N-1);
    h2 = zeros(1,N-1);
    h8 = zeros(1,N-1);
    ef = zeros(1,N-1);
    ri = zeros(1,N-1);
    sg = zeros(1,N-1);
    na = zeros(1,N-1);

    for n=1:N-1

        hn(n) = abs(sum(ev(n+1:end)));						% Hankel lower bound
        l1(n) = abs(4.0*(N+n)*sum(ev(n+1:end)));					% L1 upper bound
        h2(n) = abs(sqrt(trace(c(:,n+1:end)*diag(ev(n+1:end))*b(n+1:end,:))));	% H2 error indicator
        h8(n) = 2.0*sum(abs(ev(n+1:end)));						% Hinf upper bound
        ef(n) = 1.0 - sum(sv(1:n))/sum(sv);					% Energy error
        ri(n) = abs((4.0/pi) * atan(sqrt(sum(sv(n+1:end)))/sum(sv(1:n))));	% Robustness index
        sg(n) = abs(2.0*sum(ev(n+1:end)));						% System gain error
        na(n) = sqrt(abs(pi*sum(ev(n+1:end).^2)));				% Squareroot of area enclosed in nyquist graph
    end;

%% PLOT REDUCDED ORDER VS RELATIVE ERRORS
    if(nargin>0 && o==0), return; end; 
    figure('Name',mfilename,'NumberTitle','off');
    semilogy(1:N-1,hn,'color',[27,158,119]./255,'linewidth',2);
    hold on;
    semilogy(1:N-1,l1,'color',[217,95,2]./255,'linewidth',2);
    semilogy(1:N-1,h2,'color',[117,112,179]./255,'linewidth',2);
    semilogy(1:N-1,h8,'color',[231,41,138]./255,'linewidth',2);
    semilogy(1:N-1,sg,'color',[102,166,30]./255,'linewidth',2);
    semilogy(1:N-1,ef,'color',[230,171,2]./255,'linewidth',2);
    semilogy(1:N-1,ri,'color',[166,118,29]./255,'linewidth',2);
    semilogy(1:N-1,na,'color',[102,102,102]./255,'linewidth',2);
    hold off;
    xlim([1,N-1]);
    ylim([1e-18,1e2]);
    pbaspect([2,1.25,1]);
    legend('Hankel lower bound','L1 upper bound','H2 indicator','H8 upper bound','Energy fraction','Robustness','System gain','Nyquist area');
    set(gca,'YGrid','on');
    if(nargin>0 && o==1), print('-dsvg',[mfilename(),'.svg']); end;
end

