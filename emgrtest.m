function emgrtest(o)
% emgrtest (emgr unit tests)
% by Christian Himpe, 2015 ( http://gramian.de )
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
    J = 4;
    N = J*J;
    O = J;

    S = 0.0;
    h = 0.1;
    T = 1.0;
    L = T/h;

    A = rand(N,N); A = 0.5*(A+A'); A(1:N+1:N) = -0.55*N;
    B = rand(N,J);
    C = B';

    f = @(x,u,p) A*x + B*u + p;
    g = @(x,u,p) C*x;
    F = @(x,u,p) A'*x + C'*u;

    X = zeros(N,1);
    U = [ones(J,1),zeros(J,L-1)];
    P = rand(N,1);

    ok = @(t) fprintf([t,' - OK \n']);

    %% Default Empirical Gramians
    WC = emgr(f,g,[J,N,O],[S,h,T],'c',P); ok('Default WC');
    WO = emgr(f,g,[J,N,O],[S,h,T],'o',P); ok('Default WO');
    WX = emgr(f,g,[J,N,O],[S,h,T],'x',P); ok('Default WX');
    WY = emgr(f,F,[J,N,O],[S,h,T],'y',P); ok('Default WY');
    WS = emgr(f,g,[J,N,O],[S,h,T],'s',P); ok('Default WS');
    WI = emgr(f,g,[J,N,O],[S,h,T],'i',P); ok('Default WI');
    WJ = emgr(f,g,[J,N,O],[S,h,T],'j',P); ok('Default WJ');

    %% Parametrized Empirical Gramians
    WC = emgr(f,g,[J,N,O],[S,h,T],'c',[0.1*ones(N,1),P,ones(N,1)]); ok('Parametrized WC');
    WO = emgr(f,g,[J,N,O],[S,h,T],'o',[0.1*ones(N,1),P,ones(N,1)]); ok('Parametrized WO');
    WX = emgr(f,g,[J,N,O],[S,h,T],'x',[0.1*ones(N,1),P,ones(N,1)]); ok('Parametrized WX');
    WY = emgr(f,F,[J,N,O],[S,h,T],'y',[0.1*ones(N,1),P,ones(N,1)]); ok('Parametrized WY');
    WS = emgr(f,g,[J,N,O],[S,h,T],'s',[0.1*ones(N,1),P,ones(N,1)]); ok('Parametrized WS');
    WI = emgr(f,g,[J,N,O],[S,h,T],'i',[0.1*ones(N,1),P,ones(N,1)]); ok('Parametrized WI');
    WJ = emgr(f,g,[J,N,O],[S,h,T],'j',[0.1*ones(N,1),P,ones(N,1)]); ok('Parametrized WJ');

    %% Centering
    WC = emgr(f,g,[J,N,O],[S,h,T],'c',P,[1,0,0,0,0,0,0,0,0,0]); ok('Mean Center');
    WC = emgr(f,g,[J,N,O],[S,h,T],'c',P,[2,0,0,0,0,0,0,0,0,0]); ok('Init Center');
    WC = emgr(f,g,[J,N,O],[S,h,T],'c',P,[3,0,0,0,0,0,0,0,0,0]); ok('Steady Center');
    WC = emgr(f,g,[J,N,O],[S,h,T],'c',P,[4,0,0,0,0,0,0,0,0,0]); ok('Median Center');
    WC = emgr(f,g,[J,N,O],[S,h,T],'c',P,[5,0,0,0,0,0,0,0,0,0]); ok('Midrange Center');
    WC = emgr(f,g,[J,N,O],[S,h,T],'c',P,[6,0,0,0,0,0,0,0,0,0]); ok('RMS Center');

    %% Input Scale Spacing
    WC = emgr(f,g,[J,N,O],[S,h,T],'c',P,[0,1,0,0,0,0,0,0,0,0]); ok('Log Input Scale');
    WC = emgr(f,g,[J,N,O],[S,h,T],'c',P,[0,2,0,0,0,0,0,0,0,0]); ok('Geo Input Scale');
    WC = emgr(f,g,[J,N,O],[S,h,T],'c',P,[0,3,0,0,0,0,0,0,0,0]); ok('Single Input Scale');
    WC = emgr(f,g,[J,N,O],[S,h,T],'c',P,[0,4,0,0,0,0,0,0,0,0]); ok('Sparse Input Scale');

    %% State Scale Spacing
    WO = emgr(f,g,[J,N,O],[S,h,T],'o',P,[0,0,1,0,0,0,0,0,0,0]); ok('Log State Scale');
    WO = emgr(f,g,[J,N,O],[S,h,T],'o',P,[0,0,2,0,0,0,0,0,0,0]); ok('Geo State Scale');
    WO = emgr(f,g,[J,N,O],[S,h,T],'o',P,[0,0,3,0,0,0,0,0,0,0]); ok('Single State Scale');
    WO = emgr(f,g,[J,N,O],[S,h,T],'o',P,[0,0,4,0,0,0,0,0,0,0]); ok('Sparse State Scale');

    %% Input Rotations
    WC = emgr(f,g,[J,N,O],[S,h,T],'c',P,[0,0,0,1,0,0,0,0,0,0]); ok('Reciproce Input Rot');
    WC = emgr(f,g,[J,N,O],[S,h,T],'c',P,[0,0,0,2,0,0,0,0,0,0]); ok('Dyadic Input Rot');
    WC = emgr(f,g,[J,N,O],[S,h,T],'c',P,[0,0,0,3,0,0,0,0,0,0]); ok('Single Input Rot');

    %% State Rotations
    WO = emgr(f,g,[J,N,O],[S,h,T],'o',P,[0,0,0,0,1,0,0,0,0,0]); ok('Reciproce State Rot');
    WO = emgr(f,g,[J,N,O],[S,h,T],'o',P,[0,0,0,0,2,0,0,0,0,0]); ok('Dyadic State Rot');
    WO = emgr(f,g,[J,N,O],[S,h,T],'o',P,[0,0,0,0,3,0,0,0,0,0]); ok('Single State Rot');

    %% Double Run
    WC = emgr(f,g,[J,N,O],[S,h,T],'c',P,[0,0,0,0,0,1,0,0,0,0]); ok('Double Run WC');
    WO = emgr(f,g,[J,N,O],[S,h,T],'o',P,[0,0,0,0,0,1,0,0,0,0]); ok('Double Run WO');
    WX = emgr(f,g,[J,N,O],[S,h,T],'x',P,[0,0,0,0,0,1,0,0,0,0]); ok('Double Run WX');
    WY = emgr(f,F,[J,N,O],[S,h,T],'y',P,[0,0,0,0,0,1,0,0,0,0]); ok('Double Run WY');
    WS = emgr(f,g,[J,N,O],[S,h,T],'s',P,[0,0,0,0,0,1,0,0,0,0]); ok('Double Run WS');
    WI = emgr(f,g,[J,N,O],[S,h,T],'i',P,[0,0,0,0,0,1,0,0,0,0]); ok('Double Run WI');
    WJ = emgr(f,g,[J,N,O],[S,h,T],'j',P,[0,0,0,0,0,1,0,0,0,0]); ok('Double Run WJ');

    %% Scaled Run
    WC = emgr(f,g,[J,N,O],[S,h,T],'c',P,[0,0,0,0,0,2,0,0,0,0]); ok('Scaled Run WC');
    WO = emgr(f,g,[J,N,O],[S,h,T],'o',P,[0,0,0,0,0,2,0,0,0,0]); ok('Scaled Run WO');
    WX = emgr(f,g,[J,N,O],[S,h,T],'x',P,[0,0,0,0,0,2,0,0,0,0]); ok('Scaled Run WX');
    WY = emgr(f,F,[J,N,O],[S,h,T],'y',P,[0,0,0,0,0,2,0,0,0,0]); ok('Scaled Run WY');
    WS = emgr(f,g,[J,N,O],[S,h,T],'s',P,[0,0,0,0,0,2,0,0,0,0]); ok('Scaled Run WS');
    WI = emgr(f,g,[J,N,O],[S,h,T],'i',P,[0,0,0,0,0,2,0,0,0,0]); ok('Scaled Run WI');
    WJ = emgr(f,g,[J,N,O],[S,h,T],'j',P,[0,0,0,0,0,2,0,0,0,0]); ok('Scaled Run WJ');

    % Non-Symmetric Cross Gramians
    WX = emgr(f,g,[J,N,O],[S,h,T],'x',P,[0,0,0,0,0,0,1,0,0,0]); ok('Non-Symmetric WX');
    WJ = emgr(f,g,[J,N,O],[S,h,T],'j',P,[0,0,0,0,0,0,1,0,0,0]); ok('Non-Symmetric WJ');

    %% Robust Gramians
    WC = emgr(f,g,[J,N,O],[S,h,T],'c',P,[0,0,0,0,0,0,0,1,0,0]); ok('Robust WC');
    WY = emgr(f,F,[J,N,O],[S,h,T],'y',P,[0,0,0,0,0,0,0,1,0,0]); ok('Robust WY');

    %% Parameter Centering
    WC = emgr(f,g,[J,N,O],[S,h,T],'c',[0.1*ones(N,1),P,ones(N,1)],[0,0,0,0,0,0,0,1,1,0]); ok('ParaCent WC');
    WY = emgr(f,F,[J,N,O],[S,h,T],'y',[0.1*ones(N,1),P,ones(N,1)],[0,0,0,0,0,0,0,0,1,0]); ok('ParaCent WY');
    WY = emgr(f,g,[J,N,O],[S,h,T],'i',[0.1*ones(N,1),P,ones(N,1)],[0,0,0,0,0,0,0,0,1,0]); ok('ParaCent WI');
    WY = emgr(f,g,[J,N,O],[S,h,T],'j',[0.1*ones(N,1),P,ones(N,1)],[0,0,0,0,0,0,0,0,1,0]); ok('ParaCent WJ');

    %% Mean-Centered WS
    WS = emgr(f,g,[J,N,O],[S,h,T],'s',P,[0,0,0,0,0,0,0,0,0,1]); ok('Mean WS');

    %% Schur-Complement WI
    WI = emgr(f,g,[J,N,O],[S,h,T],'i',P,[0,0,0,0,0,0,0,0,0,1]); ok('Schur WI');

    %% Symmetric Part Gramian
    WJ = emgr(f,g,[J,N,O],[S,h,T],'x',P,[0,0,0,0,0,0,0,0,0,1]); ok('NonSym WX');

    %% Custom Solver
    ODE = @mysolver;
    WX = emgr(f,g,[J,N,O],[S,h,T],'x',P,[0,0,0,0,0,0,0,0,0,0]); ok('Custom Solver WX');
    ODE = [];

    %% Procedural Input
    WC = emgr(f,g,[J,N,O],[S,h,T],'c',P,0,@(t) exp(-0.05*t)); ok('Procedural Input WC');

    %% Chirp Input
    WC = emgr(f,g,[J,N,O],[S,h,T],'c',P,0,Inf); ok('Chirp Input WC');

    %% Demos
    disp('vernval:'); vernval

    disp('state_wx:'); state_wx
    disp('state_wy:'); state_wy
    disp('nonsym_wx'); nonsym_wx
    disp('state_bt:'); state_bt
    disp('gains_wx'); gains_wx
    disp('hierarchy:'); hierarchy
    disp('param_ws:'); param_ws
    disp('param_wi:'); param_wi
    disp('combined_wj:'); combined_wj
    disp('benchmark_ilp:'); benchmark_ilp
    disp('benchmark_lin:'); benchmark_lin
    disp('benchmark_non:'); benchmark_non
    disp('energy_wx:'); energy_wx
    disp('measure:'); measure
    disp('decentral:'); decentral
    disp('advection:'); advection
    disp('nbody:'); nbody
    disp('blackhole:'); blackhole
end
