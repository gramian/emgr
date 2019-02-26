function emgrtest()
%%% summary: emgrtest (run emgr tests)
%%% project: emgr - Empirical Gramian Framework ( https://gramian.de )
%%% authors: Christian Himpe ( 0000-0003-2194-6754 )
%%% license: BSD-2-Clause (2019)

    rand('seed',1009);								% Seed uniform random number generator

%% SYSTEM SETUP

    sys.M = 4;									% Number of inputs
    sys.N = sys.M*sys.M*sys.M;							% Number of states
    sys.Q = sys.M;								% Number of outputs
    sys.dt = 0.01;								% Time step
    sys.Tf = 1.0;								% Time horizon

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

    curios(sys,'state-reduction','linear-balanced-truncation',{[8,6,1]});
    curios(sys,'state-reduction','nonlinear-balanced-truncation',{[8,6,2]});
    curios(sys,'state-reduction','linear-direct-truncation',{[8,6,3]});
    curios(sys,'state-reduction','nonlinear-direct-truncation',{[8,6,4]});
    curios(sys,'state-reduction','linear-dominant-subspaces',{[8,6,5]});
    curios(sys,'state-reduction','nonlinear-dominant-subspaces',{[8,6,6]});

    curios(sys,'state-reduction','linear-direct-truncation',{'nonsym',[8,6,9]});
    curios(sys,'state-reduction','nonlinear-direct-truncation',{'nonsym',[8,6,10]});
    curios(sys,'state-reduction','linear-dominant-subspaces',{'nonsym',[8,6,11]});
    curios(sys,'state-reduction','nonlinear-dominant-subspaces',{'nonsym',[8,6,12]});

    curios(sys,'state-reduction','linear-balanced-truncation',{'gains',[8,6,13]});
    curios(sys,'state-reduction','linear-direct-truncation',{'gains',[8,6,15]});

    curios(sys,'state-reduction','linear-balanced-truncation',{'jacobi',[8,6,19]});
    curios(sys,'state-reduction','nonlinear-balanced-truncation',{'jacobi',[8,6,20]});
    curios(sys,'state-reduction','linear-direct-truncation',{'jacobi',[8,6,21]});
    curios(sys,'state-reduction','nonlinear-direct-truncation',{'jacobi',[8,6,22]});
    curios(sys,'state-reduction','linear-dominant-subspaces',{'jacobi',[8,6,23]});
    curios(sys,'state-reduction','nonlinear-dominant-subspaces',{'jacobi',[8,6,24]});

    curios(sys,'state-reduction','linear-balanced-truncation',{'scaled',[8,6,25]});
    curios(sys,'state-reduction','nonlinear-balanced-truncation',{'scaled',[8,6,26]});
    curios(sys,'state-reduction','linear-direct-truncation',{'scaled',[8,6,27]});
    curios(sys,'state-reduction','nonlinear-direct-truncation',{'scaled',[8,6,28]});
    curios(sys,'state-reduction','linear-dominant-subspaces',{'scaled',[8,6,29]});
    curios(sys,'state-reduction','nonlinear-dominant-subspaces',{'scaled',[8,6,30]});

    curios(sys,'state-reduction','linear-balanced-truncation',{'active',[8,6,31]});
    curios(sys,'state-reduction','nonlinear-balanced-truncation',{'active',[8,6,32]});
    curios(sys,'state-reduction','linear-direct-truncation',{'active',[8,6,33]});
    curios(sys,'state-reduction','nonlinear-direct-truncation',{'active',[8,6,34]});
    curios(sys,'state-reduction','linear-dominant-subspaces',{'active',[8,6,35]});
    curios(sys,'state-reduction','nonlinear-dominant-subspaces',{'active',[8,6,36]});

    curios(sys,'state-reduction','linear-balanced-truncation',{'tweighted',[8,6,37]});
    curios(sys,'state-reduction','nonlinear-balanced-truncation',{'tweighted',[8,6,38]});
    curios(sys,'state-reduction','linear-direct-truncation',{'tweighted',[8,6,39]});
    curios(sys,'state-reduction','nonlinear-direct-truncation',{'tweighted',[8,6,40]});
    curios(sys,'state-reduction','linear-dominant-subspaces',{'tweighted',[8,6,41]});
    curios(sys,'state-reduction','nonlinear-dominant-subspaces',{'tweighted',[8,6,42]});

    curios(sys,'state-reduction','linear-balanced-truncation',{'polynomial',[8,6,43]});
    curios(sys,'state-reduction','nonlinear-balanced-truncation',{'polynomial',[8,6,44]});
    curios(sys,'state-reduction','linear-direct-truncation',{'polynomial',[8,6,45]});
    curios(sys,'state-reduction','nonlinear-direct-truncation',{'polynomial',[8,6,46]});
    curios(sys,'state-reduction','linear-dominant-subspaces',{'polynomial',[8,6,47]});
    curios(sys,'state-reduction','nonlinear-dominant-subspaces',{'polynomial',[8,6,48]});

    % Parametric state-space reduction
    sys.p = ones(sys.N,1)*[0.01,0.5,1.0];					% Training parameter
    sys.q = rand(sys.N,5);							% Test parameter

    curios(sys,'state-reduction','linear-balanced-truncation',{[3,6,1]});
    curios(sys,'state-reduction','nonlinear-balanced-truncation',{[3,6,2]});
    curios(sys,'state-reduction','linear-direct-truncation',{[3,6,3]});
    curios(sys,'state-reduction','nonlinear-direct-truncation',{[3,6,4]});
    curios(sys,'state-reduction','linear-dominant-subspaces',{[3,6,5]});
    curios(sys,'state-reduction','nonlinear-dominant-subspaces',{[3,6,6]});

    curios(sys,'state-reduction','linear-direct-truncation',{'nonsym',[3,6,9]});
    curios(sys,'state-reduction','nonlinear-direct-truncation',{'nonsym',[3,6,10]});
    curios(sys,'state-reduction','linear-dominant-subspaces',{'nonsym',[3,6,11]});
    curios(sys,'state-reduction','nonlinear-dominant-subspaces',{'nonsym',[3,6,12]});

    curios(sys,'state-reduction','linear-direct-truncation',{'nonsym','polynomial',[3,6,15]});
    curios(sys,'state-reduction','nonlinear-direct-truncation',{'nonsym','polynomial',[3,6,16]});
    curios(sys,'state-reduction','linear-dominant-subspaces',{'nonsym','polynomial',[3,6,17]});
    curios(sys,'state-reduction','nonlinear-dominant-subspaces',{'nonsym','polynomial',[3,6,18]});

    % Parameter-space reduction tests
    sys.p = ones(sys.N,1)*[0.01,1.0];						% Training parameter range
    sys.q = rand(sys.N,5);							% Test parameter

    curios(sys,'parameter-reduction','controllability-based',{[4,4,1]});
    curios(sys,'parameter-reduction','observability-based',{[4,4,2]});
    curios(sys,'parameter-reduction','minimality-based',{[4,4,3]});
    curios(sys,'parameter-reduction','minimality-based',{'nonsym',[4,4,4]});

    curios(sys,'parameter-reduction','controllability-based',{'coarse',[4,4,5]});
    curios(sys,'parameter-reduction','observability-based',{'coarse',[4,4,6]});
    curios(sys,'parameter-reduction','minimality-based',{'coarse',[4,4,7]});
    curios(sys,'parameter-reduction','minimality-based',{'nonsym','coarse',[4,4,8]});

    curios(sys,'parameter-reduction','controllability-based',{'linpar',[4,4,9]});
    curios(sys,'parameter-reduction','observability-based',{'linpar',[4,4,10]});
    curios(sys,'parameter-reduction','minimality-based',{'linpar',[4,4,11]});
    curios(sys,'parameter-reduction','minimality-based',{'nonsym','linpar',[4,4,12]});

    curios(sys,'parameter-reduction','controllability-based',{'logpar',[4,4,13]});
    curios(sys,'parameter-reduction','observability-based',{'logpar',[4,4,14]});
    curios(sys,'parameter-reduction','minimality-based',{'logpar',[4,4,15]});
    curios(sys,'parameter-reduction','minimality-based',{'nonsym','logpar',[4,4,16]});

    % Combined state and parameter reduction
    curios(sys,'combined-reduction','controllability-based',{[1,6,1]});
    curios(sys,'combined-reduction','observability-based',{'linpar',[1,6,2]});
    curios(sys,'combined-reduction','minimality-based',{'linpar',[1,6,3]});
    curios(sys,'combined-reduction','minimality-based',{'linpar','dominant',[1,6,4]});
    curios(sys,'combined-reduction','minimality-based',{'nonsym','linpar',[1,6,5]});
    curios(sys,'combined-reduction','minimality-based',{'nonsym','linpar','dominant',[1,6,6]});

    % Distributed cross operator
    curios(sys,'state-reduction','nonlinear-direct-truncation',{[2,3,1]});
    curios(sys,'state-reduction','nonlinear-direct-truncation',{'partitioned',[2,3,4]});
    curios(sys,'parameter-reduction','minimality-based',{[2,3,2]});
    curios(sys,'parameter-reduction','minimality-based',{'partitioned',[2,3,5]});
    curios(sys,'combined-reduction','minimality-based',{[2,3,3]});
    curios(sys,'combined-reduction','minimality-based',{'partitioned',[2,3,6]});

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

    curios(sys,'decentralized-control','linear',{[1,4,1]});
    curios(sys,'decentralized-control','nonlinear',{[1,4,2]});
    curios(sys,'decentralized-control','linear',{'coherence',[1,4,3]});
    curios(sys,'decentralized-control','nonlinear',{'coherence',[1,4,4]});

% NONLINEARITY QUANTIFICATION

    % Input Nonlinearity
    sys.f = @(x,u,p,t) A*x + B*atan(u) + F*p;					% Vector field
    sys.g = @(x,u,p,t) C*x;							% Output functional

    curios(sys,'nonlinearity-quantification','input-based',{[1,4,1]});
    curios(sys,'nonlinearity-quantification','state-based',{'hold',[1,4,1]});
    curios(sys,'nonlinearity-quantification','output-based',{'hold',[1,4,1]});

    % State Nonlinearity
    sys.f = @(x,u,p,t) A*atan(x) + B*u + F*p;					% Vector field
    sys.g = @(x,u,p,t) C*x;							% Output functional

    curios(sys,'nonlinearity-quantification','input-based',{[1,4,2]});
    curios(sys,'nonlinearity-quantification','state-based',{'hold',[1,4,2]});
    curios(sys,'nonlinearity-quantification','output-based',{'hold',[1,4,2]});

    % Output Nonlinearity
    sys.f = @(x,u,p,t) A*x + B*u + F*p;						% Vector field
    sys.g = @(x,u,p,t) C*atan(x);						% Output functional

    curios(sys,'nonlinearity-quantification','input-based',{[1,4,3]});
    curios(sys,'nonlinearity-quantification','state-based',{'hold',[1,4,3]});
    curios(sys,'nonlinearity-quantification','output-based',{'hold',[1,4,3]});

    % Linear System
    sys.f = @(x,u,p,t) A*x + B*u + F*p;						% Vector field
    sys.g = @(x,u,p,t) C*x;							% Output functional

    curios(sys,'nonlinearity-quantification','input-based',{[1,4,4]});
    curios(sys,'nonlinearity-quantification','state-based',{'hold',[1,4,4]});
    curios(sys,'nonlinearity-quantification','output-based',{'hold',[1,4,4]});

%% STATE INDEX

    curios(sys,'state-index','controllability');
    curios(sys,'state-index','observability',{'hold'});
    curios(sys,'state-index','minimality',{'hold'});

%% SYSTEM INDEX

    curios(sys,'system-index','cauchy-index');
    curios(sys,'system-index','system-entropy',{'cached'});
    curios(sys,'system-index','system-gain',{'cached'});
    curios(sys,'system-index','hinf-bound',{'cached','hold'})
    curios(sys,'system-index','hankel-bound',{'cached','hold'})
    curios(sys,'system-index','nyquist-area',{'cached','hold'});
    curios(sys,'system-index','energy-fraction',{'cached','hold'});
    curios(sys,'system-index','storage-efficiency',{'cached','hold'});
    curios(sys,'system-index','ellipsoid-volume',{'cached','hold'});
    curios(sys,'system-index','io-coherence',{'cached','hold'});
end
