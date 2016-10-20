function emgrtest(m)
%%% summary: emgrtest (emgr sanity test via EOC)
%%% project: emgr - Empirical Gramian Framework ( gramian.de )
%%% authors: Christian Himpe ( 0000-0003-2194-6754 )
%%% license: 2-Clause BSD (2015--2016)
%$
    if(exist('emgr')~=2)
        error('emgr not found! Get emgr at: http://gramian.de');
    else
        global ODE;
        ODE = [];
        fprintf('emgr (version: %1.1f)\n',emgr('version'));
    end

%% Setup
    M = 4;
    N = M*M;
    Q = M;

    h1 = 0.01;
    h2 = 0.001;
    T = 1.0;

    rand('seed',1009);
    A = rand(N,N);
    A = 0.5*(A+A');
    A(1:N+1:end) = -N;
    B = rand(N,M);
    C = B';

    f = @(x,u,p,t) A*x + B*u + p;
    g = @(x,u,p,t) C*x;
    F = @(x,u,p,t) A'*x + C'*u;

    P = zeros(N,1);
    R = [zeros(N,1),ones(N,1)];

    if(nargin==0) m = ''; end

    switch(m)

        case 'oct',
            EMGR = @emgr_oct;

        case 'legacy',
            EMGR = @emgr_legacy;

        otherwise,
            EMGR = @emgr;
    end

    W = sylvester(A,A',-B*C);

    WC1 = EMGR(f,g,[M,N,Q],[h1,T],'c',P);
    WC2 = EMGR(f,g,[M,N,Q],[h2,T],'c',P);
    EOC_WC = log( norm(WC1-W,'fro')/norm(WC2-W,'fro') ) / log(h1/h2)

    WO1 = EMGR(f,g,[M,N,Q],[h1,T],'o',P);
    WO2 = EMGR(f,g,[M,N,Q],[h2,T],'o',P);
    EOC_WO = log( norm(WO1-W,'fro')/norm(WO2-W,'fro') ) / log(h1/h2)

    WX1 = EMGR(f,g,[M,N,Q],[h1,T],'x',P);
    WX2 = EMGR(f,g,[M,N,Q],[h2,T],'x',P);
    EOC_WX = log( norm(WX1-W,'fro')/norm(WX2-W,'fro') ) / log(h1/h2)

    WY1 = EMGR(f,F,[M,N,Q],[h1,T],'y',P);
    WY2 = EMGR(f,F,[M,N,Q],[h2,T],'y',P);
    EOC_WY = log( norm(WY1-W,'fro')/norm(WY2-W,'fro') ) / log(h1/h2)

    WS1 = EMGR(f,g,[M,N,Q],[h1,T],'s',R);
    WS2 = EMGR(f,g,[M,N,Q],[h2,T],'s',R);
    EOC_WS1 = log( norm(WS1{1}-W,'fro')/norm(WS2{1}-W,'fro') ) / log(h1/h2)

    WI1 = EMGR(f,g,[M,N,Q],[h1,T],'i',R);
    WI2 = EMGR(f,g,[M,N,Q],[h2,T],'i',R);
    EOC_WI1 = log( norm(WI1{1}-W,'fro')/norm(WI2{1}-W,'fro') ) / log(h1/h2)

    WJ1 = EMGR(f,g,[M,N,Q],[h1,T],'j',R);
    WJ2 = EMGR(f,g,[M,N,Q],[h2,T],'j',R);
    EOC_WJ1 = log( norm(WJ1{1}-W,'fro')/norm(WJ2{1}-W,'fro') ) / log(h1/h2)
end
