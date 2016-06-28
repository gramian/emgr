function emgrtest(o)
% emgrtest (emgr basic tests)
% by Christian Himpe, 2015-2016 ( http://gramian.de )
% released under BSD 2-Clause License ( opensource.org/licenses/BSD-2-Clause )
%*
    if(exist('emgr')~=2)
        error('emgr not found! Get emgr at: http://gramian.de');
    else
        global ODE;
        fprintf('emgr (version: %1.1f)\n',emgr('version'));
    end

%% Setup
    J = 4;
    N = J*J;
    O = J;

    h = 0.01;
    T = 1.0;
    L = floor(T/h) + 1;

    A = rand(N,N);
    A = 0.5*(A+A');
    A(1:N+1:end) = -N;
    B = rand(N,J);
    C = B';

    f = @(x,u,p) A*x + B*u + p;
    g = @(x,u,p) C*x;
    F = @(x,u,p) A'*x + C'*u;

    U = [ones(J,1),zeros(J,L-1)];
    P = zeros(N,1);
    Q = [zeros(N,1),ones(N,1)];

    if(nargin==0 || o==[] || o==0)
        WC = emgr(f,g,[J,N,O],[h,T],'c',P);
        WO = emgr(f,g,[J,N,O],[h,T],'o',P);
        WX = emgr(f,g,[J,N,O],[h,T],'x',P);
        WY = emgr(f,F,[J,N,O],[h,T],'y',P);
        WS = emgr(f,g,[J,N,O],[h,T],'s',Q);
        WI = emgr(f,g,[J,N,O],[h,T],'i',Q);
        WJ = emgr(f,g,[J,N,O],[h,T],'j',Q);
    elseif(o==1)
        WC = emgr_oct(f,g,[J,N,O],[h,T],'c',P);
        WO = emgr_oct(f,g,[J,N,O],[h,T],'o',P);
        WX = emgr_oct(f,g,[J,N,O],[h,T],'x',P);
        WY = emgr_oct(f,F,[J,N,O],[h,T],'y',P);
        WS = emgr_oct(f,g,[J,N,O],[h,T],'s',Q);
        WI = emgr_oct(f,g,[J,N,O],[h,T],'i',Q);
        WJ = emgr_oct(f,g,[J,N,O],[h,T],'j',Q);
    end

    EWCWO = norm(WC-WO,'fro')
    EWCWX = norm(WC-WX,'fro')
    EWCWY = norm(WC-WY,'fro')
    EWCWS = norm(WC-WS{1},'fro')
    EWOWI = norm(WO-WI{1},'fro')
    EWJWX = norm(WX-WJ{1},'fro')

    W = sylvester(A,A',-B*C);

    EWC = norm(WC-W,'fro')
    EWO = norm(WO-W,'fro')
    EWX = norm(WX-W,'fro')
    EWY = norm(WY-W,'fro')

    ODE = [];
end
