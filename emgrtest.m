function emgrtest(m)
%%% summary: emgrtest (run emgr tests)
%%% project: emgr - Empirical Gramian Framework ( http://gramian.de )
%%% authors: Christian Himpe ( 0000-0003-2194-6754 )
%%% license: 2-Clause BSD (2015--2018)
%$
    rand('seed',1009);								% Seed uniform random number generator

%% SYSTEM SETUP

    sys.M = 4;									% Number of inputs
    sys.N = sys.M*sys.M*sys.M;							% Number of states
    sys.Q = sys.M;								% Number of outputs
    sys.h = 0.01;								% Time step
    sys.T = 1.0;								% Time horizon

    A = -gallery('lehmer',sys.N);						% System matrix
    B = toeplitz(1:sys.N,1:sys.M);						% Input matrix
    C = B';									% Output matrix
    F = diag(10.^linspace(0,-8,sys.N));						% Source matrix

    sys.f = @(x,u,p,t) A*x + B*u + F*p;						% Vector field
    sys.g = @(x,u,p,t) C*x;							% Output functional
    sys.F = @(x,u,p,t) A'*x + C'*u;						% Adjoint vector field

%% MODEL REDUCTION

    % State-space reduction tests
    sys.p = zeros(sys.N,1);							% Training Parameters
    sys.q = zeros(sys.N,1);							% Test Parameters

    curios(sys,'state-reduction','linear-direct-truncation');
    curios(sys,'state-reduction','linear-direct-truncation',{'nonsym'});
    curios(sys,'state-reduction','nonlinear-direct-truncation');
    curios(sys,'state-reduction','nonlinear-direct-truncation',{'nonsym'});
    curios(sys,'state-reduction','linear-balanced-truncation');
    curios(sys,'state-reduction','nonlinear-balanced-truncation');

    curios(sys,'state-reduction','linear-direct-truncation',{'gains'});
    curios(sys,'state-reduction','linear-direct-truncation',{'nonsym','gains'});
    curios(sys,'state-reduction','linear-balanced-truncation',{'gains'});

    curios(sys,'state-reduction','linear-direct-truncation',{'jacobi'});
    curios(sys,'state-reduction','nonlinear-direct-truncation',{'jacobi'});
    curios(sys,'state-reduction','linear-balanced-truncation',{'jacobi'});
    curios(sys,'state-reduction','nonlinear-balanced-truncation',{'jacobi'});

    curios(sys,'state-reduction','linear-direct-truncation',{'scaled'});
    curios(sys,'state-reduction','nonlinear-direct-truncation',{'scaled'});
    curios(sys,'state-reduction','linear-balanced-truncation',{'scaled'});
    curios(sys,'state-reduction','nonlinear-balanced-truncation',{'scaled'});

    curios(sys,'state-reduction','linear-direct-truncation',{'active'});
    curios(sys,'state-reduction','nonlinear-direct-truncation',{'active'});
    curios(sys,'state-reduction','linear-balanced-truncation',{'active'});
    curios(sys,'state-reduction','nonlinear-balanced-truncation',{'active'});

    curios(sys,'state-reduction','linear-direct-truncation',{'tweighted'});
    curios(sys,'state-reduction','nonlinear-direct-truncation',{'tweighted'});
    curios(sys,'state-reduction','linear-balanced-truncation',{'tweighted'});
    curios(sys,'state-reduction','nonlinear-balanced-truncation',{'tweighted'});

    curios(sys,'state-reduction','linear-direct-truncation',{'polynomial'});
    curios(sys,'state-reduction','nonlinear-direct-truncation',{'polynomial'});
    curios(sys,'state-reduction','linear-balanced-truncation',{'polynomial'});
    curios(sys,'state-reduction','nonlinear-balanced-truncation',{'polynomial'});

    % Parametric state-space reduction
    sys.p = ones(sys.N,1)*[0.01,0.5,1.0];					% Training parameter
    sys.q = rand(sys.N,5);							% Test parameter

    curios(sys,'state-reduction','linear-direct-truncation');
    curios(sys,'state-reduction','linear-direct-truncation',{'nonsym'});
    curios(sys,'state-reduction','linear-direct-truncation',{'nonsym','polynomial'});
    curios(sys,'state-reduction','nonlinear-direct-truncation');
    curios(sys,'state-reduction','nonlinear-direct-truncation',{'nonsym'});
    curios(sys,'state-reduction','nonlinear-direct-truncation',{'nonsym','polynomial'});
    curios(sys,'state-reduction','linear-balanced-truncation');
    curios(sys,'state-reduction','nonlinear-balanced-truncation');

    % Parameter-space reduction tests
    sys.p = ones(sys.N,1)*[0.01,1.0];						% Training parameter range
    sys.q = rand(sys.N,5);							% Test parameter

    curios(sys,'parameter-reduction','controllability-based');
    curios(sys,'parameter-reduction','observability-based');
    curios(sys,'parameter-reduction','minimality-based');
    curios(sys,'parameter-reduction','minimality-based',{'nonsym'});

    curios(sys,'parameter-reduction','controllability-based',{'coarse'});
    curios(sys,'parameter-reduction','observability-based',{'coarse'});
    curios(sys,'parameter-reduction','minimality-based',{'coarse'});
    curios(sys,'parameter-reduction','minimality-based',{'nonsym','coarse'});

    curios(sys,'parameter-reduction','controllability-based',{'linpar'});
    curios(sys,'parameter-reduction','observability-based',{'linpar'});
    curios(sys,'parameter-reduction','minimality-based',{'linpar'});

    curios(sys,'parameter-reduction','controllability-based',{'logpar'});
    curios(sys,'parameter-reduction','observability-based',{'logpar'});
    curios(sys,'parameter-reduction','minimality-based',{'logpar'});

    % Combined state and parameter reduction

    curios(sys,'combined-reduction','controllability-based');
    curios(sys,'combined-reduction','observability-based',{'linpar'});
    curios(sys,'combined-reduction','minimality-based',{'linpar'});
    curios(sys,'combined-reduction','minimality-based',{'nonsym','linear'});

    % Distributed cross operator

    curios(sys,'state-reduction','nonlinear-direct-truncation',{'partitioned'});
    curios(sys,'parameter-reduction','minimality-based',{'partitioned'});
    curios(sys,'combined-reduction','minimality-based',{'partitioned'});

%% SENSITIVITY ANALYSIS

    curios(sys,'sensitivity-analysis','input-state-based');
    curios(sys,'sensitivity-analysis','input-output-based',{'hold'});

%% PARAMETER IDENTIFICATION

    curios(sys,'parameter-identification','state-output-based');
    curios(sys,'parameter-identification','input-output-based',{'hold'});
    curios(sys,'parameter-identification','state-output-based',{'coarse','hold'});
    curios(sys,'parameter-identification','input-output-based',{'coarse','hold'});

%% DECENTRALIZED CONTROL

    sys.p = zeros(sys.N,1);							% Training Parameters
    sys.q = zeros(sys.N,1);							% Test Parameters

    curios(sys,'decentralized-control','linear');
    curios(sys,'decentralized-control','nonlinear');
    curios(sys,'decentralized-control','linear',{'coherence'});
    curios(sys,'decentralized-control','nonlinear',{'coherence'});

% NONLINEARITY QUANTIFICATION

    % Input Nonlinearity
    sys.f = @(x,u,p,t) A*x + B*atan(u) + F*p;					% Vector field
    sys.g = @(x,u,p,t) C*x;							% Output functional

    curios(sys,'nonlinearity-quantification','input-based');
    curios(sys,'nonlinearity-quantification','state-based',{'hold'});
    curios(sys,'nonlinearity-quantification','output-based',{'hold'});

    % State Nonlinearity
    sys.f = @(x,u,p,t) A*atan(x) + B*u + F*p;					% Vector field
    sys.g = @(x,u,p,t) C*x;							% Output functional

    curios(sys,'nonlinearity-quantification','input-based');
    curios(sys,'nonlinearity-quantification','state-based',{'hold'});
    curios(sys,'nonlinearity-quantification','output-based',{'hold'});

    % Output Nonlinearity
    sys.f = @(x,u,p,t) A*x + B*u + F*p;						% Vector field
    sys.g = @(x,u,p,t) C*atan(x);						% Output functional

    curios(sys,'nonlinearity-quantification','input-based');
    curios(sys,'nonlinearity-quantification','state-based',{'hold'});
    curios(sys,'nonlinearity-quantification','output-based',{'hold'});

    % Linear System
    sys.f = @(x,u,p,t) A*x + B*u + F*p;						% Vector field
    sys.g = @(x,u,p,t) C*x;							% Output functional

    curios(sys,'nonlinearity-quantification','input-based');
    curios(sys,'nonlinearity-quantification','state-based',{'hold'});
    curios(sys,'nonlinearity-quantification','output-based',{'hold'});

%% STATE INDEX

    curios(sys,'state-index','controllability');
    curios(sys,'state-index','observability',{'hold'});
    curios(sys,'state-index','minimality',{'hold'});

%% SYSTEM INDEX

    curios(sys,'system-index','cauchy-index');
    curios(sys,'system-index','system-entropy');
    curios(sys,'system-index','system-gain');
    curios(sys,'system-index','hinf-bound',{'hold'})
    curios(sys,'system-index','hankel-bound',{'hold'})
    curios(sys,'system-index','system-symmetry',{'hold'});
    curios(sys,'system-index','nyquist-area',{'hold'});
    curios(sys,'system-index','storage-efficiency',{'hold'});
    curios(sys,'system-index','robustness-index',{'hold'});
end
