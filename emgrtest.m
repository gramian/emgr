function emgrtest(o)
% emgrtest (emgr unit tests)
% by Christian Himpe, 2015 ( http://gramian.de )
% released under BSD 2-Clause License ( opensource.org/licenses/BSD-2-Clause )
%*
    if(exist('emgr')~=2)
        error('emgr not found! Get emgr at: http://gramian.de');
    else
        global ODE; ODE = [];
        fprintf('emgr (version: %g)\n',emgr('version'));
    end

    %% Setup
    J = 4;
    N = J*J;
    O = J;

    h = 0.1;
    T = 1.0;
    L = floor(T/h) + 1;

    A = rand(N,N); A = 0.5*(A+A'); A(1:N+1:N) = -0.55*N;
    B = rand(N,J);
    C = B';

    f = @(x,u,p) A*x + B*u + p;
    g = @(x,u,p) C*x;
    F = @(x,u,p) A'*x + C'*u;

    X = zeros(N,1);
    U = [ones(J,1),zeros(J,L-1)];
    P = 0.5*rand(N,1)+0.25;
    Q = ones(N,1)*[0.25,0.75];

    ok = @(t) fprintf([t,' - OK \n']);

    %% Default Empirical Gramians
    WC = emgr(f,g,[J,N,O],[h,T],'c',P); ok('Default WC');
    WO = emgr(f,g,[J,N,O],[h,T],'o',P); ok('Default WO');
    WX = emgr(f,g,[J,N,O],[h,T],'x',P); ok('Default WX');
    WY = emgr(f,F,[J,N,O],[h,T],'y',P); ok('Default WY');
    WS = emgr(f,g,[J,N,O],[h,T],'s',Q); ok('Default WS');
    WI = emgr(f,g,[J,N,O],[h,T],'i',Q); ok('Default WI');
    WJ = emgr(f,g,[J,N,O],[h,T],'j',Q); ok('Default WJ');

    %% Parametrized Empirical Gramians
    WC = emgr(f,g,[J,N,O],[h,T],'c',[zeros(N,1),ones(N,1)]); ok('Parametrized WC');
    WO = emgr(f,g,[J,N,O],[h,T],'o',[zeros(N,1),ones(N,1)]); ok('Parametrized WO');
    WX = emgr(f,g,[J,N,O],[h,T],'x',[zeros(N,1),ones(N,1)]); ok('Parametrized WX');
    WY = emgr(f,F,[J,N,O],[h,T],'y',[zeros(N,1),ones(N,1)]); ok('Parametrized WY');

    %% Centering
    WC = emgr(f,g,[J,N,O],[h,T],'c',P,[1,0,0,0,0,0,0,0,0,0,0,0]); ok('Mean Center');
    WC = emgr(f,g,[J,N,O],[h,T],'c',P,[2,0,0,0,0,0,0,0,0,0,0,0]); ok('Init Center');
    WC = emgr(f,g,[J,N,O],[h,T],'c',P,[3,0,0,0,0,0,0,0,0,0,0,0]); ok('Steady Center');
    WC = emgr(f,g,[J,N,O],[h,T],'c',P,[4,0,0,0,0,0,0,0,0,0,0,0]); ok('Median Center');
    WC = emgr(f,g,[J,N,O],[h,T],'c',P,[5,0,0,0,0,0,0,0,0,0,0,0]); ok('Midrange Center');
    WC = emgr(f,g,[J,N,O],[h,T],'c',P,[6,0,0,0,0,0,0,0,0,0,0,0]); ok('RMS Center');

    %% Input Scale Spacing
    WC = emgr(f,g,[J,N,O],[h,T],'c',P,[0,1,0,0,0,0,0,0,0,0,0,0]); ok('Log Input Scale');
    WC = emgr(f,g,[J,N,O],[h,T],'c',P,[0,2,0,0,0,0,0,0,0,0,0,0]); ok('Geo Input Scale');
    WC = emgr(f,g,[J,N,O],[h,T],'c',P,[0,3,0,0,0,0,0,0,0,0,0,0]); ok('Single Input Scale');
    WC = emgr(f,g,[J,N,O],[h,T],'c',P,[0,4,0,0,0,0,0,0,0,0,0,0]); ok('Sparse Input Scale');

    %% State Scale Spacing
    WO = emgr(f,g,[J,N,O],[h,T],'o',P,[0,0,1,0,0,0,0,0,0,0,0,0]); ok('Log State Scale');
    WO = emgr(f,g,[J,N,O],[h,T],'o',P,[0,0,2,0,0,0,0,0,0,0,0,0]); ok('Geo State Scale');
    WO = emgr(f,g,[J,N,O],[h,T],'o',P,[0,0,3,0,0,0,0,0,0,0,0,0]); ok('Single State Scale');
    WO = emgr(f,g,[J,N,O],[h,T],'o',P,[0,0,4,0,0,0,0,0,0,0,0,0]); ok('Sparse State Scale');

    %% Input Rotations
    WC = emgr(f,g,[J,N,O],[h,T],'c',P,[0,0,0,1,0,0,0,0,0,0,0,0]); ok('Reciproce Input Rot');
    WC = emgr(f,g,[J,N,O],[h,T],'c',P,[0,0,0,2,0,0,0,0,0,0,0,0]); ok('Dyadic Input Rot');
    WC = emgr(f,g,[J,N,O],[h,T],'c',P,[0,0,0,3,0,0,0,0,0,0,0,0]); ok('Single Input Rot');

    %% State Rotations
    WO = emgr(f,g,[J,N,O],[h,T],'o',P,[0,0,0,0,1,0,0,0,0,0,0,0]); ok('Reciproce State Rot');
    WO = emgr(f,g,[J,N,O],[h,T],'o',P,[0,0,0,0,2,0,0,0,0,0,0,0]); ok('Dyadic State Rot');
    WO = emgr(f,g,[J,N,O],[h,T],'o',P,[0,0,0,0,3,0,0,0,0,0,0,0]); ok('Single State Rot');

    %% Double Run
    WC = emgr(f,g,[J,N,O],[h,T],'c',P,[0,0,0,0,0,1,0,0,0,0,0,0]); ok('Double Run WC');
    WO = emgr(f,g,[J,N,O],[h,T],'o',P,[0,0,0,0,0,1,0,0,0,0,0,0]); ok('Double Run WO');
    WX = emgr(f,g,[J,N,O],[h,T],'x',P,[0,0,0,0,0,1,0,0,0,0,0,0]); ok('Double Run WX');
    WY = emgr(f,F,[J,N,O],[h,T],'y',P,[0,0,0,0,0,1,0,0,0,0,0,0]); ok('Double Run WY');
    WS = emgr(f,g,[J,N,O],[h,T],'s',Q,[0,0,0,0,0,1,0,0,0,0,0,0]); ok('Double Run WS');
    WI = emgr(f,g,[J,N,O],[h,T],'i',Q,[0,0,0,0,0,1,0,0,0,0,0,0]); ok('Double Run WI');
    WJ = emgr(f,g,[J,N,O],[h,T],'j',Q,[0,0,0,0,0,1,0,0,0,0,0,0]); ok('Double Run WJ');

    %% Scaled Run
    WC = emgr(f,g,[J,N,O],[h,T],'c',P,[0,0,0,0,0,2,0,0,0,0,0,0]); ok('Scaled Run WC');
    WO = emgr(f,g,[J,N,O],[h,T],'o',P,[0,0,0,0,0,2,0,0,0,0,0,0]); ok('Scaled Run WO');
    WX = emgr(f,g,[J,N,O],[h,T],'x',P,[0,0,0,0,0,2,0,0,0,0,0,0]); ok('Scaled Run WX');
    WY = emgr(f,F,[J,N,O],[h,T],'y',P,[0,0,0,0,0,2,0,0,0,0,0,0]); ok('Scaled Run WY');
    WS = emgr(f,g,[J,N,O],[h,T],'s',Q,[0,0,0,0,0,2,0,0,0,0,0,0]); ok('Scaled Run WS');
    WI = emgr(f,g,[J,N,O],[h,T],'i',Q,[0,0,0,0,0,2,0,0,0,0,0,0]); ok('Scaled Run WI');
    WJ = emgr(f,g,[J,N,O],[h,T],'j',Q,[0,0,0,0,0,2,0,0,0,0,0,0]); ok('Scaled Run WJ');

    % Non-Symmetric Cross Gramians
    WX = emgr(f,g,[J,N,O],[h,T],'x',P,[0,0,0,0,0,0,1,0,0,0,0,0]); ok('Non-Symmetric WX');
    WJ = emgr(f,g,[J,N,O],[h,T],'j',Q,[0,0,0,0,0,0,1,0,0,0,0,0]); ok('Non-Symmetric WJ');

    %% Robust Gramians
    WC = emgr(f,g,[J,N,O],[h,T],'c',Q,[0,0,0,0,0,0,0,1,0,0,0,0]); ok('Robust WC');
    WY = emgr(f,F,[J,N,O],[h,T],'y',Q,[0,0,0,0,0,0,0,1,0,0,0,0]); ok('Robust WY');

    % Active Parameters


    %% Linear Parameter Centering
    WY = emgr(f,F,[J,N,O],[h,T],'s',[zeros(N,1),ones(N,1)],[0,0,0,0,0,0,0,0,1,0,0,0]); ok('Lin ParaCent WS');
    WY = emgr(f,g,[J,N,O],[h,T],'i',[zeros(N,1),ones(N,1)],[0,0,0,0,0,0,0,0,1,0,0,0]); ok('Lin ParaCent WI');
    WY = emgr(f,g,[J,N,O],[h,T],'j',[zeros(N,1),ones(N,1)],[0,0,0,0,0,0,0,0,1,0,0,0]); ok('Lin ParaCent WJ');

    %% Logarithmic Parameter Centering
    WY = emgr(f,F,[J,N,O],[h,T],'s',[zeros(N,1),ones(N,1)],[0,0,0,0,0,0,0,0,2,0,0,0]); ok('Log ParaCent WS');
    WY = emgr(f,g,[J,N,O],[h,T],'i',[zeros(N,1),ones(N,1)],[0,0,0,0,0,0,0,0,2,0,0,0]); ok('Log ParaCent WI');
    WY = emgr(f,g,[J,N,O],[h,T],'j',[zeros(N,1),ones(N,1)],[0,0,0,0,0,0,0,0,2,0,0,0]); ok('Log ParaCent WJ');

    %% Mean-Centered WS
    WS = emgr(f,g,[J,N,O],[h,T],'s',Q,[0,0,0,0,0,0,0,0,0,1,0,0]); ok('Mean WS');

    %% Schur-Complement WI
    WI = emgr(f,g,[J,N,O],[h,T],'i',Q,[0,0,0,0,0,0,0,0,0,1,0,0]); ok('Schur WI');

    %% Symmetric Part Gramian
    WJ = emgr(f,g,[J,N,O],[h,T],'x',Q,[0,0,0,0,0,0,0,0,0,1,0,0]); ok('Detailed Schur WX');

    %% Custom Solver
    ODE = @mysolver;
    WX = emgr(f,g,[J,N,O],[h,T],'x',P,[0,0,0,0,0,0,0,0,0,0,0,0]); ok('Custom Solver WX');
    ODE = [];

    %% Procedural Input
    WC = emgr(f,g,[J,N,O],[h,T],'c',P,0,@(t) exp(-0.05*t)); ok('Procedural Input WC');

    %% Chirp Input
    WC = emgr(f,g,[J,N,O],[h,T],'c',P,0,Inf); ok('Chirp Input WC');

    %% Demos
    disp('vernval:'); vernval

    disp('state_wx:'); state_wx(0)
    disp('state_wy:'); state_wy(0)
    disp('nonsym_wx'); nonsym_wx(0)
    disp('state_bt:'); state_bt(0)
    disp('gains_wx'); gains_wx(0)
    disp('hierarchy:'); hierarchy(0)
    disp('param_ws:'); param_ws(0)
    disp('param_wi:'); param_wi(0)
    disp('combined_wj:'); combined_wj(0)
    disp('benchmark_ilp:'); benchmark_ilp(0)
    disp('benchmark_lin:'); benchmark_lin(0)
    disp('benchmark_non:'); benchmark_non(0)
    disp('energy_wx:'); energy_wx(0)
    disp('measure:'); measure(0)
    disp('decentral:'); decentral(0)
    disp('advection:'); advection(0)
    disp('nbody:'); nbody(0)
    disp('blackhole:'); blackhole(0)
end
