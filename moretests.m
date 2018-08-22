function moretests()
%%% summary: moretests
%%% project: emgr - EMpirical GRamian Framework ( https://gramian.de )
%%% authors: Christian Himpe ( 0000-0003-2194-6754 )
%%% license: BSD-2-Clause (2018)
%$
    if(exist('emgr')~=2)
        error('emgr not found! Get emgr at: http://gramian.de');
    else
        global ODE;
        ODE = [];
        fprintf('emgr (version: %1.1f)\n',emgr('version'));
    end

%% SYSTEM SETUP
    M = 4;				% number of inputs
    N = M*M;				% number of states
    Q = M;				% number of outputs
    h = 0.01;				% time step size
    T = 1.0;				% time horizon
    X = linspace(0,1,N)';		% initial state
    U = @(t) ones(M,1)*(t<=h)/h;	% impulse input function
    P = 0.5+0.5*cos(1:N)';		% parameter
    R = [0.5*ones(N,1),ones(N,1)];	% parameter range

    A = diag(N*(sin(1:N)-1));		% system matrix
    B = sin(toeplitz(1:N,1:M)) + 1.0;	% input matrix
    C = cos(toeplitz(1:M,1:N)) + 1.0;	% output matrix

    LIN = @(x,u,p,t) A*x + B*u + p;	% vector field
    ADJ = @(x,u,p,t) A'*x + C'*u;	% adjoint vector field
    OUT = @(x,u,p,t) C*x;		% output functional

    % figure; plot(0:h:T,ODE(LIN,OUT,[h,T],X,U,P)); return;

    if(exist('OCTAVE_VERSION','builtin'))
        EMGR = @emgr_oct;
    elseif(verLessThan('matlab','9.1'))
        EMGR = @emgr_lgc;
    else
        EMGR = @emgr;
    end

%% Controllability Gramian

    fprintf('Empirical Controllability Gramian: ');
    W = {}; k = 1;

    % Baseline
    W{k} = svd(EMGR(LIN,OUT,[M,N,Q],[h,T],'c',P,[0,0,0,0,0,0,0,0,0,0,0,0],1,0,X)); fprintf('#'); k = k + 1;

    % Centering
    W{k} = svd(EMGR(LIN,OUT,[M,N,Q],[h,T],'c',P,[1,0,0,0,0,0,0,0,0,0,0,0],1,0,X)); fprintf('#'); k = k + 1;
    W{k} = svd(EMGR(LIN,OUT,[M,N,Q],[h,T],'c',P,[2,0,0,0,0,0,0,0,0,0,0,0],1,0,X)); fprintf('#'); k = k + 1;
    W{k} = svd(EMGR(LIN,OUT,[M,N,Q],[h,T],'c',P,[3,0,0,0,0,0,0,0,0,0,0,0],1,0,X)); fprintf('#'); k = k + 1;
    W{k} = svd(EMGR(LIN,OUT,[M,N,Q],[h,T],'c',P,[4,0,0,0,0,0,0,0,0,0,0,0],1,0,X)); fprintf('#'); k = k + 1;
    W{k} = svd(EMGR(LIN,OUT,[M,N,Q],[h,T],'c',P,[5,0,0,0,0,0,0,0,0,0,0,0],1,0,X)); fprintf('#'); k = k + 1;

    % Input scales
    W{k} = svd(EMGR(LIN,OUT,[M,N,Q],[h,T],'c',P,[0,1,0,0,0,0,0,0,0,0,0,0],1,0,X)); fprintf('#'); k = k + 1;
    W{k} = svd(EMGR(LIN,OUT,[M,N,Q],[h,T],'c',P,[0,2,0,0,0,0,0,0,0,0,0,0],1,0,X)); fprintf('#'); k = k + 1;
    W{k} = svd(EMGR(LIN,OUT,[M,N,Q],[h,T],'c',P,[0,3,0,0,0,0,0,0,0,0,0,0],1,0,X)); fprintf('#'); k = k + 1;
    W{k} = svd(EMGR(LIN,OUT,[M,N,Q],[h,T],'c',P,[0,4,0,0,0,0,0,0,0,0,0,0],1,0,X)); fprintf('#'); k = k + 1;

    % Input Rotations
    W{k} = svd(EMGR(LIN,OUT,[M,N,Q],[h,T],'c',P,[0,0,0,1,0,0,0,0,0,0,0,0],1,0,X)); fprintf('#'); k = k + 1;

    % Normalization
    W{k} = svd(EMGR(LIN,OUT,[M,N,Q],[h,T],'c',P,[0,0,0,0,0,1,0,0,0,0,0,0],1,0,X)); fprintf('#'); k = k + 1;
    W{k} = svd(EMGR(LIN,OUT,[M,N,Q],[h,T],'c',P,[0,0,0,0,0,2,0,0,0,0,0,0],1,0,X)); fprintf('#'); k = k + 1;

    plotsingvals(W,[2,4,1],'Controllability Gramian'); fprintf('\n');

%% Observability Gramian

    fprintf('Empirical Observability Gramian: ');
    W = {}; k = 1;

    % Baseline
    W{k} = svd(EMGR(LIN,OUT,[M,N,Q],[h,T],'o',P,[0,0,0,0,0,0,0,0,0,0,0,0],1,0,X)); fprintf('#'); k = k + 1;

    % Centering
    W{k} = svd(EMGR(LIN,OUT,[M,N,Q],[h,T],'o',P,[1,0,0,0,0,0,0,0,0,0,0,0],1,0,X)); fprintf('#'); k = k + 1;
    W{k} = svd(EMGR(LIN,OUT,[M,N,Q],[h,T],'o',P,[2,0,0,0,0,0,0,0,0,0,0,0],1,0,X)); fprintf('#'); k = k + 1;
    W{k} = svd(EMGR(LIN,OUT,[M,N,Q],[h,T],'o',P,[3,0,0,0,0,0,0,0,0,0,0,0],1,0,X)); fprintf('#'); k = k + 1;
    W{k} = svd(EMGR(LIN,OUT,[M,N,Q],[h,T],'o',P,[4,0,0,0,0,0,0,0,0,0,0,0],1,0,X)); fprintf('#'); k = k + 1;
    W{k} = svd(EMGR(LIN,OUT,[M,N,Q],[h,T],'o',P,[5,0,0,0,0,0,0,0,0,0,0,0],1,0,X)); fprintf('#'); k = k + 1;

    % State scales
    W{k} = svd(EMGR(LIN,OUT,[M,N,Q],[h,T],'o',P,[0,0,1,0,0,0,0,0,0,0,0,0],1,0,X)); fprintf('#'); k = k + 1;
    W{k} = svd(EMGR(LIN,OUT,[M,N,Q],[h,T],'o',P,[0,0,2,0,0,0,0,0,0,0,0,0],1,0,X)); fprintf('#'); k = k + 1;
    W{k} = svd(EMGR(LIN,OUT,[M,N,Q],[h,T],'o',P,[0,0,3,0,0,0,0,0,0,0,0,0],1,0,X)); fprintf('#'); k = k + 1;
    W{k} = svd(EMGR(LIN,OUT,[M,N,Q],[h,T],'o',P,[0,0,4,0,0,0,0,0,0,0,0,0],1,0,X)); fprintf('#'); k = k + 1;

    % State Rotations
    W{k} = svd(EMGR(LIN,OUT,[M,N,Q],[h,T],'o',P,[0,0,0,0,1,0,0,0,0,0,0,0],1,0,X)); fprintf('#'); k = k + 1;

    % Normalization
    W{k} = svd(EMGR(LIN,OUT,[M,N,Q],[h,T],'o',P,[0,0,0,0,0,1,0,0,0,0,0,0],1,0,X)); fprintf('#'); k = k + 1;
    W{k} = svd(EMGR(LIN,OUT,[M,N,Q],[h,T],'o',P,[0,0,0,0,0,1,0,0,0,0,0,0],1,0,X)); fprintf('#'); k = k + 1;

    % Extra Input
    W{k} = svd(EMGR(LIN,OUT,[M,N,Q],[h,T],'o',P,[0,0,0,0,0,0,0,1,0,0,0,0],1,0,X)); fprintf('#'); k = k + 1;

    plotsingvals(W,[2,4,2],'Observability Gramian'); fprintf('\n');

%% Cross Gramian

    fprintf('Empirical Cross Gramian: ');
    W = {}; k = 1;

    % Baseline
    W{k} = svd(EMGR(LIN,OUT,[M,N,Q],[h,T],'x',P,[0,0,0,0,0,0,0,0,0,0,0,0],1,0,X)); fprintf('#'); k = k + 1;

    % Centering
    W{k} = svd(EMGR(LIN,OUT,[M,N,Q],[h,T],'x',P,[1,0,0,0,0,0,0,0,0,0,0,0],1,0,X)); fprintf('#'); k = k + 1;
    W{k} = svd(EMGR(LIN,OUT,[M,N,Q],[h,T],'x',P,[2,0,0,0,0,0,0,0,0,0,0,0],1,0,X)); fprintf('#'); k = k + 1;
    W{k} = svd(EMGR(LIN,OUT,[M,N,Q],[h,T],'x',P,[3,0,0,0,0,0,0,0,0,0,0,0],1,0,X)); fprintf('#'); k = k + 1;
    W{k} = svd(EMGR(LIN,OUT,[M,N,Q],[h,T],'x',P,[4,0,0,0,0,0,0,0,0,0,0,0],1,0,X)); fprintf('#'); k = k + 1;
    W{k} = svd(EMGR(LIN,OUT,[M,N,Q],[h,T],'x',P,[5,0,0,0,0,0,0,0,0,0,0,0],1,0,X)); fprintf('#'); k = k + 1;

    % Input scales
    W{k} = svd(EMGR(LIN,OUT,[M,N,Q],[h,T],'x',P,[0,1,0,0,0,0,0,0,0,0,0,0],1,0,X)); fprintf('#'); k = k + 1;
    W{k} = svd(EMGR(LIN,OUT,[M,N,Q],[h,T],'x',P,[0,2,0,0,0,0,0,0,0,0,0,0],1,0,X)); fprintf('#'); k = k + 1;
    W{k} = svd(EMGR(LIN,OUT,[M,N,Q],[h,T],'x',P,[0,3,0,0,0,0,0,0,0,0,0,0],1,0,X)); fprintf('#'); k = k + 1;
    W{k} = svd(EMGR(LIN,OUT,[M,N,Q],[h,T],'x',P,[0,4,0,0,0,0,0,0,0,0,0,0],1,0,X)); fprintf('#'); k = k + 1;

    % State scales
    W{k} = svd(EMGR(LIN,OUT,[M,N,Q],[h,T],'o',P,[0,0,1,0,0,0,0,0,0,0,0,0],1,0,X)); fprintf('#'); k = k + 1;
    W{k} = svd(EMGR(LIN,OUT,[M,N,Q],[h,T],'o',P,[0,0,2,0,0,0,0,0,0,0,0,0],1,0,X)); fprintf('#'); k = k + 1;
    W{k} = svd(EMGR(LIN,OUT,[M,N,Q],[h,T],'o',P,[0,0,3,0,0,0,0,0,0,0,0,0],1,0,X)); fprintf('#'); k = k + 1;
    W{k} = svd(EMGR(LIN,OUT,[M,N,Q],[h,T],'o',P,[0,0,4,0,0,0,0,0,0,0,0,0],1,0,X)); fprintf('#'); k = k + 1;

    % Input Rotations
    W{k} = svd(EMGR(LIN,OUT,[M,N,Q],[h,T],'x',P,[0,0,0,1,0,0,0,0,0,0,0,0],1,0,X)); fprintf('#'); k = k + 1;

    % State Rotations
    W{k} = svd(EMGR(LIN,OUT,[M,N,Q],[h,T],'x',P,[0,0,0,0,1,0,0,0,0,0,0,0],1,0,X)); fprintf('#'); k = k + 1;

    % Normalization
    W{k} = svd(EMGR(LIN,OUT,[M,N,Q],[h,T],'x',P,[0,0,0,0,0,1,0,0,0,0,0,0],1,0,X)); fprintf('#'); k = k + 1;
    W{k} = svd(EMGR(LIN,OUT,[M,N,Q],[h,T],'x',P,[0,0,0,0,0,1,0,0,0,0,0,0],1,0,X)); fprintf('#'); k = k + 1;

    % Non-Symmetric Cross Gramian
    W{k} = svd(EMGR(LIN,OUT,[M,N,Q],[h,T],'x',P,[0,0,0,0,0,0,1,0,0,0,0,0],1,0,X)); fprintf('#'); k = k + 1; % fast drop

    % Extra Input
    W{k} = svd(EMGR(LIN,OUT,[M,N,Q],[h,T],'x',P,[0,0,0,0,0,0,0,1,0,0,0,0],1,0,X)); fprintf('#'); k = k + 1;

    plotsingvals(W,[2,4,3],'Cross Gramian'); fprintf('\n');

%% Linear Cross Gramian

    fprintf('Empirical Linear Cross Gramian: ');
    W = {}; k = 1;

    % Baseline
    W{k} = svd(EMGR(LIN,ADJ,[M,N,Q],[h,T],'y',P,[0,0,0,0,0,0,0,0,0,0,0,0],1,0,X)); fprintf('#'); k = k + 1;

    % Centering
    W{k} = svd(EMGR(LIN,ADJ,[M,N,Q],[h,T],'y',P,[1,0,0,0,0,0,0,0,0,0,0,0],1,0,X)); fprintf('#'); k = k + 1;
    W{k} = svd(EMGR(LIN,ADJ,[M,N,Q],[h,T],'y',P,[2,0,0,0,0,0,0,0,0,0,0,0],1,0,X)); fprintf('#'); k = k + 1;
    W{k} = svd(EMGR(LIN,ADJ,[M,N,Q],[h,T],'y',P,[3,0,0,0,0,0,0,0,0,0,0,0],1,0,X)); fprintf('#'); k = k + 1;
    W{k} = svd(EMGR(LIN,ADJ,[M,N,Q],[h,T],'y',P,[4,0,0,0,0,0,0,0,0,0,0,0],1,0,X)); fprintf('#'); k = k + 1;
    W{k} = svd(EMGR(LIN,ADJ,[M,N,Q],[h,T],'y',P,[5,0,0,0,0,0,0,0,0,0,0,0],1,0,X)); fprintf('#'); k = k + 1;

    % Input scales
    W{k} = svd(EMGR(LIN,ADJ,[M,N,Q],[h,T],'y',P,[0,1,0,0,0,0,0,0,0,0,0,0],1,0,X)); fprintf('#'); k = k + 1;
    W{k} = svd(EMGR(LIN,ADJ,[M,N,Q],[h,T],'y',P,[0,2,0,0,0,0,0,0,0,0,0,0],1,0,X)); fprintf('#'); k = k + 1;
    W{k} = svd(EMGR(LIN,ADJ,[M,N,Q],[h,T],'y',P,[0,3,0,0,0,0,0,0,0,0,0,0],1,0,X)); fprintf('#'); k = k + 1;
    W{k} = svd(EMGR(LIN,ADJ,[M,N,Q],[h,T],'y',P,[0,4,0,0,0,0,0,0,0,0,0,0],1,0,X)); fprintf('#'); k = k + 1;

    % Input Rotations
    W{k} = svd(EMGR(LIN,ADJ,[M,N,Q],[h,T],'y',P,[0,0,0,1,0,0,0,0,0,0,0,0],1,0,X)); fprintf('#'); k = k + 1;

    % Normalization
    W{k} = svd(EMGR(LIN,ADJ,[M,N,Q],[h,T],'y',P,[0,0,0,0,0,1,0,0,0,0,0,0],1,0,X)); fprintf('#'); k = k + 1;
    W{k} = svd(EMGR(LIN,ADJ,[M,N,Q],[h,T],'y',P,[0,0,0,0,0,1,0,0,0,0,0,0],1,0,X)); fprintf('#'); k = k + 1;

    % Non-Symmetric Cross Gramian
    W{k} = svd(EMGR(LIN,ADJ,[M,N,Q],[h,T],'y',P,[0,0,0,0,0,0,1,0,0,0,0,0],1,0,X)); fprintf('#'); k = k + 1; % fast drop

    plotsingvals(W,[2,4,4],'Linear Cross Gramian'); fprintf('\n');

%% Sensitivity Gramian

    fprintf('Empirical Sensitivity Gramian: ');
    W = {}; k = 1;

    % Baseline
    w = EMGR(LIN,OUT,[M,N,Q],[h,T],'s',R,[0,0,0,0,0,0,0,0,0,0,0,0],1,0,X); W{k} = sort(w{2},'descend'); fprintf('#'); k = k + 1;

    % Centering
    w = EMGR(LIN,OUT,[M,N,Q],[h,T],'s',R,[1,0,0,0,0,0,0,0,0,0,0,0],1,0,X); W{k} = sort(w{2},'descend'); fprintf('#'); k = k + 1;
    w = EMGR(LIN,OUT,[M,N,Q],[h,T],'s',R,[2,0,0,0,0,0,0,0,0,0,0,0],1,0,X); W{k} = sort(w{2},'descend'); fprintf('#'); k = k + 1;
    w = EMGR(LIN,OUT,[M,N,Q],[h,T],'s',R,[3,0,0,0,0,0,0,0,0,0,0,0],1,0,X); W{k} = sort(w{2},'descend'); fprintf('#'); k = k + 1;
    w = EMGR(LIN,OUT,[M,N,Q],[h,T],'s',R,[4,0,0,0,0,0,0,0,0,0,0,0],1,0,X); W{k} = sort(w{2},'descend'); fprintf('#'); k = k + 1;
    w = EMGR(LIN,OUT,[M,N,Q],[h,T],'s',R,[5,0,0,0,0,0,0,0,0,0,0,0],1,0,X); W{k} = sort(w{2},'descend'); fprintf('#'); k = k + 1;

    % Input scales
    w = EMGR(LIN,OUT,[M,N,Q],[h,T],'s',R,[0,1,0,0,0,0,0,0,0,0,0,0],1,0,X); W{k} = sort(w{2},'descend'); fprintf('#'); k = k + 1;
    w = EMGR(LIN,OUT,[M,N,Q],[h,T],'s',R,[0,2,0,0,0,0,0,0,0,0,0,0],1,0,X); W{k} = sort(w{2},'descend'); fprintf('#'); k = k + 1;
    w = EMGR(LIN,OUT,[M,N,Q],[h,T],'s',R,[0,3,0,0,0,0,0,0,0,0,0,0],1,0,X); W{k} = sort(w{2},'descend'); fprintf('#'); k = k + 1;
    w = EMGR(LIN,OUT,[M,N,Q],[h,T],'s',R,[0,4,0,0,0,0,0,0,0,0,0,0],1,0,X); W{k} = sort(w{2},'descend'); fprintf('#'); k = k + 1;

    % Input Rotations
    w = EMGR(LIN,OUT,[M,N,Q],[h,T],'s',R,[0,0,0,1,0,0,0,0,0,0,0,0],1,0,X); W{k} = sort(w{2},'descend'); fprintf('#'); k = k + 1;

    % Extra Input
    w = EMGR(LIN,OUT,[M,N,Q],[h,T],'s',R,[0,0,0,0,0,0,0,1,0,0,0,0],1,0,X); W{k} = sort(w{2},'descend'); fprintf('#'); k = k + 1;

    % Parameter Scaling
    w = EMGR(LIN,OUT,[M,N,Q],[h,T],'s',R,[0,0,0,0,0,0,0,0,1,0,0,0],1,0,X); W{k} = sort(w{2},'descend'); fprintf('#'); k = k + 1;
    w = EMGR(LIN,OUT,[M,N,Q],[h,T],'s',R,[0,0,0,0,0,0,0,0,2,0,0,0],1,0,X); W{k} = sort(w{2},'descend'); fprintf('#'); k = k + 1;

    plotsingvals(W,[2,4,5],'Sensitivity Gramian');
    W = {}; k = 1;

    % Averaging Type
    w = EMGR(LIN,OUT,[M,N,Q],[h,T],'s',R,[0,0,0,0,0,0,0,0,0,1,0,0],1,0,X); W{k} = sort(w{2},'descend'); fprintf('#'); k = k + 1;

    % Centering
    w = EMGR(LIN,OUT,[M,N,Q],[h,T],'s',R,[1,0,0,0,0,0,0,0,0,1,0,0],1,0,X); W{k} = sort(w{2},'descend'); fprintf('#'); k = k + 1; % fast drop
    w = EMGR(LIN,OUT,[M,N,Q],[h,T],'s',R,[2,0,0,0,0,0,0,0,0,1,0,0],1,0,X); W{k} = sort(w{2},'descend'); fprintf('#'); k = k + 1;
    w = EMGR(LIN,OUT,[M,N,Q],[h,T],'s',R,[3,0,0,0,0,0,0,0,0,1,0,0],1,0,X); W{k} = sort(w{2},'descend'); fprintf('#'); k = k + 1;
    w = EMGR(LIN,OUT,[M,N,Q],[h,T],'s',R,[4,0,0,0,0,0,0,0,0,1,0,0],1,0,X); W{k} = sort(w{2},'descend'); fprintf('#'); k = k + 1;
    w = EMGR(LIN,OUT,[M,N,Q],[h,T],'s',R,[5,0,0,0,0,0,0,0,0,1,0,0],1,0,X); W{k} = sort(w{2},'descend'); fprintf('#'); k = k + 1;

    % Input scales
    w = EMGR(LIN,OUT,[M,N,Q],[h,T],'s',R,[0,1,0,0,0,0,0,0,0,1,0,0],1,0,X); W{k} = sort(w{2},'descend'); fprintf('#'); k = k + 1;
    w = EMGR(LIN,OUT,[M,N,Q],[h,T],'s',R,[0,2,0,0,0,0,0,0,0,1,0,0],1,0,X); W{k} = sort(w{2},'descend'); fprintf('#'); k = k + 1;
    w = EMGR(LIN,OUT,[M,N,Q],[h,T],'s',R,[0,3,0,0,0,0,0,0,0,1,0,0],1,0,X); W{k} = sort(w{2},'descend'); fprintf('#'); k = k + 1;
    w = EMGR(LIN,OUT,[M,N,Q],[h,T],'s',R,[0,4,0,0,0,0,0,0,0,1,0,0],1,0,X); W{k} = sort(w{2},'descend'); fprintf('#'); k = k + 1;

    % State scales
    w = EMGR(LIN,OUT,[M,N,Q],[h,T],'s',R,[0,0,1,0,0,0,0,0,0,1,0,0],1,0,X); W{k} = sort(w{2},'descend'); fprintf('#'); k = k + 1;
    w = EMGR(LIN,OUT,[M,N,Q],[h,T],'s',R,[0,0,2,0,0,0,0,0,0,1,0,0],1,0,X); W{k} = sort(w{2},'descend'); fprintf('#'); k = k + 1;
    w = EMGR(LIN,OUT,[M,N,Q],[h,T],'s',R,[0,0,3,0,0,0,0,0,0,1,0,0],1,0,X); W{k} = sort(w{2},'descend'); fprintf('#'); k = k + 1;
    w = EMGR(LIN,OUT,[M,N,Q],[h,T],'s',R,[0,0,4,0,0,0,0,0,0,1,0,0],1,0,X); W{k} = sort(w{2},'descend'); fprintf('#'); k = k + 1;

    % Input Rotations
    w = EMGR(LIN,OUT,[M,N,Q],[h,T],'s',R,[0,0,0,1,0,0,0,0,0,1,0,0],1,0,X); W{k} = sort(w{2},'descend'); fprintf('#'); k = k + 1;

    % State Rotations
    w = EMGR(LIN,OUT,[M,N,Q],[h,T],'s',R,[0,0,0,0,1,0,0,0,0,1,0,0],1,0,X); W{k} = sort(w{2},'descend'); fprintf('#'); k = k + 1;

    % Extra Input
    w = EMGR(LIN,OUT,[M,N,Q],[h,T],'s',R,[0,0,0,0,0,0,0,1,0,1,0,0],1,0,X); W{k} = sort(w{2},'descend'); fprintf('#'); k = k + 1;

    % Parameter Scaling
    w = EMGR(LIN,OUT,[M,N,Q],[h,T],'s',R,[0,0,0,0,0,0,0,0,1,1,0,0],1,0,X); W{k} = sort(w{2},'descend'); fprintf('#'); k = k + 1;
    w = EMGR(LIN,OUT,[M,N,Q],[h,T],'s',R,[0,0,0,0,0,0,0,0,2,1,0,0],1,0,X); W{k} = sort(w{2},'descend'); fprintf('#'); k = k + 1;

    plotsingvals(W,[2,4,6],'Sensitivity Gramian (Input-Output)'); fprintf('\n');

%% Identifiability Gramian

    fprintf('Empirical Identifiability Gramian: ');
    W = {}; k = 1;

    % Baseline
    w = EMGR(LIN,OUT,[M,N,Q],[h,T],'i',R,[0,0,0,0,0,0,0,0,0,0,0,0]); W{k} = svd(w{2}); fprintf('#'); k = k + 1;

    % Centering
    w = EMGR(LIN,OUT,[M,N,Q],[h,T],'i',R,[1,0,0,0,0,0,0,0,0,0,0,0]); W{k} = svd(w{2}); fprintf('#'); k = k + 1;
    w = EMGR(LIN,OUT,[M,N,Q],[h,T],'i',R,[2,0,0,0,0,0,0,0,0,0,0,0]); W{k} = svd(w{2}); fprintf('#'); k = k + 1;
    w = EMGR(LIN,OUT,[M,N,Q],[h,T],'i',R,[3,0,0,0,0,0,0,0,0,0,0,0]); W{k} = svd(w{2}); fprintf('#'); k = k + 1;
    w = EMGR(LIN,OUT,[M,N,Q],[h,T],'i',R,[4,0,0,0,0,0,0,0,0,0,0,0]); W{k} = svd(w{2}); fprintf('#'); k = k + 1;
    w = EMGR(LIN,OUT,[M,N,Q],[h,T],'i',R,[5,0,0,0,0,0,0,0,0,0,0,0]); W{k} = svd(w{2}); fprintf('#'); k = k + 1;

    % State scales
    w = EMGR(LIN,OUT,[M,N,Q],[h,T],'i',R,[0,0,1,0,0,0,0,0,0,0,0,0]); W{k} = svd(w{2}); fprintf('#'); k = k + 1;
    w = EMGR(LIN,OUT,[M,N,Q],[h,T],'i',R,[0,0,2,0,0,0,0,0,0,0,0,0]); W{k} = svd(w{2}); fprintf('#'); k = k + 1;
    w = EMGR(LIN,OUT,[M,N,Q],[h,T],'i',R,[0,0,3,0,0,0,0,0,0,0,0,0]); W{k} = svd(w{2}); fprintf('#'); k = k + 1;
    w = EMGR(LIN,OUT,[M,N,Q],[h,T],'i',R,[0,0,4,0,0,0,0,0,0,0,0,0]); W{k} = svd(w{2}); fprintf('#'); k = k + 1;

    % State Rotations
    w = EMGR(LIN,OUT,[M,N,Q],[h,T],'i',R,[0,0,0,0,1,0,0,0,0,0,0,0]); W{k} = svd(w{2}); fprintf('#'); k = k + 1;

    % Extra Input
    w = EMGR(LIN,OUT,[M,N,Q],[h,T],'i',R,[0,0,0,0,0,0,0,1,0,0,0,0]); W{k} = svd(w{2}); fprintf('#'); k = k + 1;

    % Parameter Scaling
    w = EMGR(LIN,OUT,[M,N,Q],[h,T],'i',R,[0,0,0,0,0,0,0,0,1,0,0,0]); W{k} = svd(w{2}); fprintf('#'); k = k + 1;
    w = EMGR(LIN,OUT,[M,N,Q],[h,T],'i',R,[0,0,0,0,0,0,0,0,2,0,0,0]); W{k} = svd(w{2}); fprintf('#'); k = k + 1;

    % Schur Complement
    w = EMGR(LIN,OUT,[M,N,Q],[h,T],'i',R,[0,0,0,0,0,0,0,0,0,1,0,0]); W{k} = svd(w{2}); fprintf('#'); k = k + 1;

    plotsingvals(W,[2,4,7],'Identifiability Gramian'); fprintf('\n');

%% Joint Gramian

    fprintf('Empirical Joint Gramian: ');
    W = {}; k = 1;

    % Baseline
    w = EMGR(LIN,OUT,[M,N,Q],[h,T],'j',R,[0,0,0,0,0,0,0,0,0,0,0,0]); W{k} = svd(w{2}); fprintf('#'); k = k + 1;

    % Centering
    w = EMGR(LIN,OUT,[M,N,Q],[h,T],'j',R,[1,0,0,0,0,0,0,0,0,0,0,0]); W{k} = svd(w{2}); fprintf('#'); k = k + 1;
    w = EMGR(LIN,OUT,[M,N,Q],[h,T],'j',R,[2,0,0,0,0,0,0,0,0,0,0,0]); W{k} = svd(w{2}); fprintf('#'); k = k + 1;
    w = EMGR(LIN,OUT,[M,N,Q],[h,T],'j',R,[3,0,0,0,0,0,0,0,0,0,0,0]); W{k} = svd(w{2}); fprintf('#'); k = k + 1;
    w = EMGR(LIN,OUT,[M,N,Q],[h,T],'j',R,[4,0,0,0,0,0,0,0,0,0,0,0]); W{k} = svd(w{2}); fprintf('#'); k = k + 1;
    w = EMGR(LIN,OUT,[M,N,Q],[h,T],'j',R,[5,0,0,0,0,0,0,0,0,0,0,0]); W{k} = svd(w{2}); fprintf('#'); k = k + 1;

    % Input scales
    w = EMGR(LIN,OUT,[M,N,Q],[h,T],'j',R,[0,1,0,0,0,0,0,0,0,0,0,0]); W{k} = svd(w{2}); fprintf('#'); k = k + 1;
    w = EMGR(LIN,OUT,[M,N,Q],[h,T],'j',R,[0,2,0,0,0,0,0,0,0,0,0,0]); W{k} = svd(w{2}); fprintf('#'); k = k + 1;
    w = EMGR(LIN,OUT,[M,N,Q],[h,T],'j',R,[0,3,0,0,0,0,0,0,0,0,0,0]); W{k} = svd(w{2}); fprintf('#'); k = k + 1;
    w = EMGR(LIN,OUT,[M,N,Q],[h,T],'j',R,[0,4,0,0,0,0,0,0,0,0,0,0]); W{k} = svd(w{2}); fprintf('#'); k = k + 1;

    % State scales
    w = EMGR(LIN,OUT,[M,N,Q],[h,T],'j',R,[0,0,1,0,0,0,0,0,0,0,0,0]); W{k} = svd(w{2}); fprintf('#'); k = k + 1;
    w = EMGR(LIN,OUT,[M,N,Q],[h,T],'j',R,[0,0,2,0,0,0,0,0,0,0,0,0]); W{k} = svd(w{2}); fprintf('#'); k = k + 1;
    w = EMGR(LIN,OUT,[M,N,Q],[h,T],'j',R,[0,0,3,0,0,0,0,0,0,0,0,0]); W{k} = svd(w{2}); fprintf('#'); k = k + 1;
    w = EMGR(LIN,OUT,[M,N,Q],[h,T],'j',R,[0,0,4,0,0,0,0,0,0,0,0,0]); W{k} = svd(w{2}); fprintf('#'); k = k + 1;

    % Input Rotations
    w = EMGR(LIN,OUT,[M,N,Q],[h,T],'j',R,[0,0,0,1,0,0,0,0,0,0,0,0]); W{k} = svd(w{2}); fprintf('#'); k = k + 1;

    % State Rotations
    w = EMGR(LIN,OUT,[M,N,Q],[h,T],'j',R,[0,0,0,0,1,0,0,0,0,0,0,0]); W{k} = svd(w{2}); fprintf('#'); k = k + 1;

    % Non-Symmetric Cross Gramian
    w = EMGR(LIN,OUT,[M,N,Q],[h,T],'j',R,[0,0,0,0,0,0,1,0,0,0,0,0]); W{k} = svd(w{2}); fprintf('#'); k = k + 1; % fast drop

    % Extra Input
    w = EMGR(LIN,OUT,[M,N,Q],[h,T],'j',R,[0,0,0,0,0,0,0,1,0,0,0,0]); W{k} = svd(w{2}); fprintf('#'); k = k + 1;

    % Parameter Scaling
    w = EMGR(LIN,OUT,[M,N,Q],[h,T],'j',R,[0,0,0,0,0,0,0,0,1,0,0,0]); W{k} = svd(w{2}); fprintf('#'); k = k + 1;
    w = EMGR(LIN,OUT,[M,N,Q],[h,T],'j',R,[0,0,0,0,0,0,0,0,2,0,0,0]); W{k} = svd(w{2}); fprintf('#'); k = k + 1;

    % Schur Complement
    w = EMGR(LIN,OUT,[M,N,Q],[h,T],'j',R,[0,0,0,0,0,0,0,0,0,1,0,0]); W{k} = svd(w{2}); fprintf('#'); k = k + 1;

    plotsingvals(W,[2,4,8],'Cross-Identifiability Gramian'); fprintf('\n\n');
end

function plotsingvals(W,s,l)

    if(s(3)==1),figure; end;
    subplot(s(1),s(2),s(3))
    c = winter(numel(W)-1);
    semilogy(W{2}./max(W{2}),'Linewidth',2,'color',c(1,:));
    hold on;
    for k=3:numel(W)
        semilogy(W{k}./max(W{k}),'Linewidth',2,'color',c(k-1,:));
    end
    semilogy(W{1}./max(W{1}),'Linewidth',2,'color',[1,0,0]);
    hold off;
    xlim([1,numel(W{1})]);
    xlabel(l)
    legend(strsplit(num2str([2:numel(W),1]))','location','eastoutside');
end
