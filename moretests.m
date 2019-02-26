function moretests(ut)
%%% summary: moretests
%%% project: emgr - EMpirical GRamian Framework ( https://gramian.de )
%%% authors: Christian Himpe ( 0000-0003-2194-6754 )
%%% license: BSD-2-Clause (2019)

    if( (nargin<1)  || isempty(ut) ), ut = 1; end

    if(exist('OCTAVE_VERSION','builtin'))
        EMGR = @emgr_oct;
        fprintf('emgr (version: %1.1f%s)\n',EMGR('version'),'-oct');
    elseif(verLessThan('matlab','9.1'))
        EMGR = @emgr_lgc;
        fprintf('emgr (version: %1.1f%s)\n',EMGR('version'),'-lgc');
    else
        EMGR = @emgr;
        fprintf('emgr (version: %1.1f%s)\n',EMGR('version'),'');
    end

%% SYSTEM SETUP
    M = 4;				% number of inputs
    N = M*M;				% number of states
    Q = M;				% number of outputs
    dt = 0.01;				% time step size
    Tf = 1.0;				% time horizon
    X0 = linspace(0,1,N)';		% initial state
    P = 0.5+0.5*cos(1:N)';		% parameter
    R = [0.5*ones(N,1),ones(N,1)];	% parameter range

    A = diag(N*(sin(1:N)-1));		% system matrix
    B = sin(toeplitz(1:N,1:M)) + 1.0;	% input matrix
    C = cos(toeplitz(1:M,1:N)) + 1.0;	% output matrix

    LIN = @(x,u,p,t) A*x + B*u + p;	% vector field
    ADJ = @(x,u,p,t) A'*x + C'*u;	% adjoint vector field
    OUT = @(x,u,p,t) C*x;		% output functional

%% Controllability Gramian

    fprintf('Empirical Controllability Gramian: ');
    W = {};

    % Baseline
    W{end+1} = svd(EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'c',P,[0,0,0,0,0,0,0,0,0,0,0,0],ut,0,X0)); fprintf('#');

    % Centering
    W{end+1} = svd(EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'c',P,[1,0,0,0,0,0,0,0,0,0,0,0],ut,0,X0)); fprintf('#');
    W{end+1} = svd(EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'c',P,[2,0,0,0,0,0,0,0,0,0,0,0],ut,0,X0)); fprintf('#');
    W{end+1} = svd(EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'c',P,[3,0,0,0,0,0,0,0,0,0,0,0],ut,0,X0)); fprintf('#');
    W{end+1} = svd(EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'c',P,[4,0,0,0,0,0,0,0,0,0,0,0],ut,0,X0)); fprintf('#');
    W{end+1} = svd(EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'c',P,[5,0,0,0,0,0,0,0,0,0,0,0],ut,0,X0)); fprintf('#');
    W{end+1} = svd(EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'c',P,[6,0,0,0,0,0,0,0,0,0,0,0],ut,0,X0)); fprintf('#');

    % Input scales
    W{end+1} = svd(EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'c',P,[0,1,0,0,0,0,0,0,0,0,0,0],ut,0,X0)); fprintf('#');
    W{end+1} = svd(EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'c',P,[0,2,0,0,0,0,0,0,0,0,0,0],ut,0,X0)); fprintf('#');
    W{end+1} = svd(EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'c',P,[0,3,0,0,0,0,0,0,0,0,0,0],ut,0,X0)); fprintf('#'); % fast drop
    W{end+1} = svd(EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'c',P,[0,4,0,0,0,0,0,0,0,0,0,0],ut,0,X0)); fprintf('#');

    % Input Rotations
    W{end+1} = svd(EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'c',P,[0,0,0,1,0,0,0,0,0,0,0,0],ut,0,X0)); fprintf('#');

    % Normalization
    W{end+1} = svd(EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'c',P,[0,0,0,0,0,1,0,0,0,0,0,0],ut,0,X0)); fprintf('#');
    W{end+1} = svd(EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'c',P,[0,0,0,0,0,2,0,0,0,0,0,0],ut,0,X0)); fprintf('#');

    plotsingvals(W,[2,6,1],'Controllability Gramian',10); fprintf('\n');

%% Observability Gramian

    fprintf('Empirical Observability Gramian:   ');
    W = {};

    % Baseline
    W{end+1} = svd(EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'o',P,[0,0,0,0,0,0,0,0,0,0,0,0],ut,0,X0)); fprintf('#');

    % Centering
    W{end+1} = svd(EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'o',P,[1,0,0,0,0,0,0,0,0,0,0,0],ut,0,X0)); fprintf('#');
    W{end+1} = svd(EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'o',P,[2,0,0,0,0,0,0,0,0,0,0,0],ut,0,X0)); fprintf('#');
    W{end+1} = svd(EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'o',P,[3,0,0,0,0,0,0,0,0,0,0,0],ut,0,X0)); fprintf('#');
    W{end+1} = svd(EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'o',P,[4,0,0,0,0,0,0,0,0,0,0,0],ut,0,X0)); fprintf('#');
    W{end+1} = svd(EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'o',P,[5,0,0,0,0,0,0,0,0,0,0,0],ut,0,X0)); fprintf('#');
    W{end+1} = svd(EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'o',P,[6,0,0,0,0,0,0,0,0,0,0,0],ut,0,X0)); fprintf('#');

    % State scales
    W{end+1} = svd(EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'o',P,[0,0,1,0,0,0,0,0,0,0,0,0],ut,0,X0)); fprintf('#');
    W{end+1} = svd(EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'o',P,[0,0,2,0,0,0,0,0,0,0,0,0],ut,0,X0)); fprintf('#');
    W{end+1} = svd(EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'o',P,[0,0,3,0,0,0,0,0,0,0,0,0],ut,0,X0)); fprintf('#'); % fast drop
    W{end+1} = svd(EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'o',P,[0,0,4,0,0,0,0,0,0,0,0,0],ut,0,X0)); fprintf('#');

    % State Rotations
    W{end+1} = svd(EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'o',P,[0,0,0,0,1,0,0,0,0,0,0,0],ut,0,X0)); fprintf('#');

    % Normalization
    W{end+1} = svd(EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'o',P,[0,0,0,0,0,1,0,0,0,0,0,0],ut,0,X0)); fprintf('#');
    W{end+1} = svd(EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'o',P,[0,0,0,0,0,2,0,0,0,0,0,0],ut,0,X0)); fprintf('#');

    % Average Observability Gramian
    W{end+1} = svd(EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'o',P,[0,0,0,0,0,0,1,0,0,0,0,0],ut,0,X0)); fprintf('#');

    % Extra Input
    W{end+1} = svd(EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'o',P,[0,0,0,0,0,0,0,1,0,0,0,0],ut,0,X0)); fprintf('#');

    plotsingvals(W,[2,6,2],'Observability Gramian',10); fprintf('\n');

%% Linear Cross Gramian

    fprintf('Empirical Linear Cross Gramian:    ');
    W = {};

    % Baseline
    W{end+1} = svd(EMGR(LIN,ADJ,[M,N,Q],[dt,Tf],'y',P,[0,0,0,0,0,0,0,0,0,0,0,0],ut,0,X0)); fprintf('#');

    % Centering
    W{end+1} = svd(EMGR(LIN,ADJ,[M,N,Q],[dt,Tf],'y',P,[1,0,0,0,0,0,0,0,0,0,0,0],ut,0,X0)); fprintf('#');
    W{end+1} = svd(EMGR(LIN,ADJ,[M,N,Q],[dt,Tf],'y',P,[2,0,0,0,0,0,0,0,0,0,0,0],ut,0,X0)); fprintf('#');
    W{end+1} = svd(EMGR(LIN,ADJ,[M,N,Q],[dt,Tf],'y',P,[3,0,0,0,0,0,0,0,0,0,0,0],ut,0,X0)); fprintf('#');
    W{end+1} = svd(EMGR(LIN,ADJ,[M,N,Q],[dt,Tf],'y',P,[4,0,0,0,0,0,0,0,0,0,0,0],ut,0,X0)); fprintf('#');
    W{end+1} = svd(EMGR(LIN,ADJ,[M,N,Q],[dt,Tf],'y',P,[5,0,0,0,0,0,0,0,0,0,0,0],ut,0,X0)); fprintf('#');
    W{end+1} = svd(EMGR(LIN,ADJ,[M,N,Q],[dt,Tf],'y',P,[6,0,0,0,0,0,0,0,0,0,0,0],ut,0,X0)); fprintf('#');

    % Scales
    W{end+1} = svd(EMGR(LIN,ADJ,[M,N,Q],[dt,Tf],'y',P,[0,1,0,0,0,0,0,0,0,0,0,0],ut,0,X0)); fprintf('#');
    W{end+1} = svd(EMGR(LIN,ADJ,[M,N,Q],[dt,Tf],'y',P,[0,2,0,0,0,0,0,0,0,0,0,0],ut,0,X0)); fprintf('#');
    W{end+1} = svd(EMGR(LIN,ADJ,[M,N,Q],[dt,Tf],'y',P,[0,3,0,0,0,0,0,0,0,0,0,0],ut,0,X0)); fprintf('#'); % fast drop
    W{end+1} = svd(EMGR(LIN,ADJ,[M,N,Q],[dt,Tf],'y',P,[0,4,0,0,0,0,0,0,0,0,0,0],ut,0,X0)); fprintf('#');

    % Rotations
    W{end+1} = svd(EMGR(LIN,ADJ,[M,N,Q],[dt,Tf],'y',P,[0,0,0,1,0,0,0,0,0,0,0,0],ut,0,X0)); fprintf('#');

    % Normalization
    W{end+1} = svd(EMGR(LIN,ADJ,[M,N,Q],[dt,Tf],'y',P,[0,0,0,0,0,1,0,0,0,0,0,0],ut,0,X0)); fprintf('#');
    W{end+1} = svd(EMGR(LIN,ADJ,[M,N,Q],[dt,Tf],'y',P,[0,0,0,0,0,2,0,0,0,0,0,0],ut,0,X0)); fprintf('#');

    plotsingvals(W,[2,6,3],'Linear Cross Gramian',10);

    W = {};

    % Non-Symmetric Cross Gramian
    W{end+1} = svd(EMGR(LIN,ADJ,[M,N,Q],[dt,Tf],'y',P,[0,0,0,0,0,0,1,0,0,0,0,0],ut,0,X0)); fprintf('#');

    % Centering
    W{end+1} = svd(EMGR(LIN,ADJ,[M,N,Q],[dt,Tf],'y',P,[1,0,0,0,0,0,1,0,0,0,0,0],ut,0,X0)); fprintf('#');
    W{end+1} = svd(EMGR(LIN,ADJ,[M,N,Q],[dt,Tf],'y',P,[2,0,0,0,0,0,1,0,0,0,0,0],ut,0,X0)); fprintf('#');
    W{end+1} = svd(EMGR(LIN,ADJ,[M,N,Q],[dt,Tf],'y',P,[3,0,0,0,0,0,1,0,0,0,0,0],ut,0,X0)); fprintf('#');
    W{end+1} = svd(EMGR(LIN,ADJ,[M,N,Q],[dt,Tf],'y',P,[4,0,0,0,0,0,1,0,0,0,0,0],ut,0,X0)); fprintf('#');
    W{end+1} = svd(EMGR(LIN,ADJ,[M,N,Q],[dt,Tf],'y',P,[5,0,0,0,0,0,1,0,0,0,0,0],ut,0,X0)); fprintf('#');
    W{end+1} = svd(EMGR(LIN,ADJ,[M,N,Q],[dt,Tf],'y',P,[6,0,0,0,0,0,1,0,0,0,0,0],ut,0,X0)); fprintf('#');

    % Scales
    W{end+1} = svd(EMGR(LIN,ADJ,[M,N,Q],[dt,Tf],'y',P,[0,1,0,0,0,0,1,0,0,0,0,0],ut,0,X0)); fprintf('#');
    W{end+1} = svd(EMGR(LIN,ADJ,[M,N,Q],[dt,Tf],'y',P,[0,2,0,0,0,0,1,0,0,0,0,0],ut,0,X0)); fprintf('#');
    W{end+1} = svd(EMGR(LIN,ADJ,[M,N,Q],[dt,Tf],'y',P,[0,3,0,0,0,0,1,0,0,0,0,0],ut,0,X0)); fprintf('#'); % fast drop
    W{end+1} = svd(EMGR(LIN,ADJ,[M,N,Q],[dt,Tf],'y',P,[0,4,0,0,0,0,1,0,0,0,0,0],ut,0,X0)); fprintf('#');

    % Rotations
    W{end+1} = svd(EMGR(LIN,ADJ,[M,N,Q],[dt,Tf],'y',P,[0,0,0,1,0,0,1,0,0,0,0,0],ut,0,X0)); fprintf('#'); % deep drop

    % Normalization
    W{end+1} = svd(EMGR(LIN,ADJ,[M,N,Q],[dt,Tf],'y',P,[0,0,0,0,0,1,1,0,0,0,0,0],ut,0,X0)); fprintf('#');
    W{end+1} = svd(EMGR(LIN,ADJ,[M,N,Q],[dt,Tf],'y',P,[0,0,0,0,0,2,1,0,0,0,0,0],ut,0,X0)); fprintf('#');

    plotsingvals(W,[2,6,4],'Linear Cross Gramian (Non-Symmetric)',12); fprintf('\n');

%% Cross Gramian

    fprintf('Empirical Cross Gramian:           ');
    W = {};

    % Baseline
    W{end+1} = svd(EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'x',P,[0,0,0,0,0,0,0,0,0,0,0,0],ut,0,X0)); fprintf('#');

    % Centering
    W{end+1} = svd(EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'x',P,[1,0,0,0,0,0,0,0,0,0,0,0],ut,0,X0)); fprintf('#');
    W{end+1} = svd(EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'x',P,[2,0,0,0,0,0,0,0,0,0,0,0],ut,0,X0)); fprintf('#');
    W{end+1} = svd(EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'x',P,[3,0,0,0,0,0,0,0,0,0,0,0],ut,0,X0)); fprintf('#');
    W{end+1} = svd(EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'x',P,[4,0,0,0,0,0,0,0,0,0,0,0],ut,0,X0)); fprintf('#');
    W{end+1} = svd(EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'x',P,[5,0,0,0,0,0,0,0,0,0,0,0],ut,0,X0)); fprintf('#');
    W{end+1} = svd(EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'x',P,[6,0,0,0,0,0,0,0,0,0,0,0],ut,0,X0)); fprintf('#');

    % Input scales
    W{end+1} = svd(EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'x',P,[0,1,0,0,0,0,0,0,0,0,0,0],ut,0,X0)); fprintf('#');
    W{end+1} = svd(EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'x',P,[0,2,0,0,0,0,0,0,0,0,0,0],ut,0,X0)); fprintf('#');
    W{end+1} = svd(EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'x',P,[0,3,0,0,0,0,0,0,0,0,0,0],ut,0,X0)); fprintf('#');
    W{end+1} = svd(EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'x',P,[0,4,0,0,0,0,0,0,0,0,0,0],ut,0,X0)); fprintf('#');

    % State scales
    W{end+1} = svd(EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'x',P,[0,0,1,0,0,0,0,0,0,0,0,0],ut,0,X0)); fprintf('#');
    W{end+1} = svd(EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'x',P,[0,0,2,0,0,0,0,0,0,0,0,0],ut,0,X0)); fprintf('#');
    W{end+1} = svd(EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'x',P,[0,0,3,0,0,0,0,0,0,0,0,0],ut,0,X0)); fprintf('#');
    W{end+1} = svd(EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'x',P,[0,0,4,0,0,0,0,0,0,0,0,0],ut,0,X0)); fprintf('#');

    % Input Rotations
    W{end+1} = svd(EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'x',P,[0,0,0,1,0,0,0,0,0,0,0,0],ut,0,X0)); fprintf('#');

    % State Rotations
    W{end+1} = svd(EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'x',P,[0,0,0,0,1,0,0,0,0,0,0,0],ut,0,X0)); fprintf('#'); % fast drop

    % Normalization
    W{end+1} = svd(EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'x',P,[0,0,0,0,0,1,0,0,0,0,0,0],ut,0,X0)); fprintf('#');
    W{end+1} = svd(EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'x',P,[0,0,0,0,0,2,0,0,0,0,0,0],ut,0,X0)); fprintf('#');

    % Extra Input
    W{end+1} = svd(EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'x',P,[0,0,0,0,0,0,0,1,0,0,0,0],ut,0,X0)); fprintf('#');

    plotsingvals(W,[2,6,5],'Cross Gramian',17);
    W = {};

    % Baseline
    W{end+1} = svd(EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'x',P,[0,0,0,0,0,0,1,0,0,0,0,0],ut,0,X0)); fprintf('#');

    % Centering
    W{end+1} = svd(EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'x',P,[1,0,0,0,0,0,1,0,0,0,0,0],ut,0,X0)); fprintf('#');
    W{end+1} = svd(EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'x',P,[2,0,0,0,0,0,1,0,0,0,0,0],ut,0,X0)); fprintf('#');
    W{end+1} = svd(EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'x',P,[3,0,0,0,0,0,1,0,0,0,0,0],ut,0,X0)); fprintf('#');
    W{end+1} = svd(EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'x',P,[4,0,0,0,0,0,1,0,0,0,0,0],ut,0,X0)); fprintf('#');
    W{end+1} = svd(EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'x',P,[5,0,0,0,0,0,1,0,0,0,0,0],ut,0,X0)); fprintf('#');
    W{end+1} = svd(EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'x',P,[6,0,0,0,0,0,1,0,0,0,0,0],ut,0,X0)); fprintf('#');

    % Input scales
    W{end+1} = svd(EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'x',P,[0,1,0,0,0,0,1,0,0,0,0,0],ut,0,X0)); fprintf('#');
    W{end+1} = svd(EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'x',P,[0,2,0,0,0,0,1,0,0,0,0,0],ut,0,X0)); fprintf('#');
    W{end+1} = svd(EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'x',P,[0,3,0,0,0,0,1,0,0,0,0,0],ut,0,X0)); fprintf('#');
    W{end+1} = svd(EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'x',P,[0,4,0,0,0,0,1,0,0,0,0,0],ut,0,X0)); fprintf('#');

    % State scales
    W{end+1} = svd(EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'x',P,[0,0,1,0,0,0,1,0,0,0,0,0],ut,0,X0)); fprintf('#');
    W{end+1} = svd(EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'x',P,[0,0,2,0,0,0,1,0,0,0,0,0],ut,0,X0)); fprintf('#');
    W{end+1} = svd(EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'x',P,[0,0,3,0,0,0,1,0,0,0,0,0],ut,0,X0)); fprintf('#');
    W{end+1} = svd(EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'x',P,[0,0,4,0,0,0,1,0,0,0,0,0],ut,0,X0)); fprintf('#');

    % Input Rotations
    W{end+1} = svd(EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'x',P,[0,0,0,1,0,0,1,0,0,0,0,0],ut,0,X0)); fprintf('#');

    % State Rotations
    W{end+1} = svd(EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'x',P,[0,0,0,0,1,0,1,0,0,0,0,0],ut,0,X0)); fprintf('#'); % fast drop

    % Normalization
    W{end+1} = svd(EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'x',P,[0,0,0,0,0,1,1,0,0,0,0,0],ut,0,X0)); fprintf('#');
    W{end+1} = svd(EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'x',P,[0,0,0,0,0,2,1,0,0,0,0,0],ut,0,X0)); fprintf('#'); 

    % Extra Input
    W{end+1} = svd(EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'x',P,[0,0,0,0,0,0,1,1,0,0,0,0],ut,0,X0)); fprintf('#');

    plotsingvals(W,[2,6,6],'Cross Gramian (Non-Symmetric)',17); fprintf('\n');

%% Sensitivity Gramian

    fprintf('Empirical Sensitivity Gramian:     ');
    W = {};

    % Baseline
    w = EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'s',R,[0,0,0,0,0,0,0,0,0,0,0,0],ut,0,X0); W{end+1} = sort(w{2},'descend'); fprintf('#');

    % Centering
    w = EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'s',R,[1,0,0,0,0,0,0,0,0,0,0,0],ut,0,X0); W{end+1} = sort(w{2},'descend'); fprintf('#');
    w = EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'s',R,[2,0,0,0,0,0,0,0,0,0,0,0],ut,0,X0); W{end+1} = sort(w{2},'descend'); fprintf('#');
    w = EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'s',R,[3,0,0,0,0,0,0,0,0,0,0,0],ut,0,X0); W{end+1} = sort(w{2},'descend'); fprintf('#');
    w = EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'s',R,[4,0,0,0,0,0,0,0,0,0,0,0],ut,0,X0); W{end+1} = sort(w{2},'descend'); fprintf('#');
    w = EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'s',R,[5,0,0,0,0,0,0,0,0,0,0,0],ut,0,X0); W{end+1} = sort(w{2},'descend'); fprintf('#');
    w = EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'s',R,[6,0,0,0,0,0,0,0,0,0,0,0],ut,0,X0); W{end+1} = sort(w{2},'descend'); fprintf('#');

    % Input scales
    w = EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'s',R,[0,1,0,0,0,0,0,0,0,0,0,0],ut,0,X0); W{end+1} = sort(w{2},'descend'); fprintf('#');
    w = EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'s',R,[0,2,0,0,0,0,0,0,0,0,0,0],ut,0,X0); W{end+1} = sort(w{2},'descend'); fprintf('#');
    w = EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'s',R,[0,3,0,0,0,0,0,0,0,0,0,0],ut,0,X0); W{end+1} = sort(w{2},'descend'); fprintf('#');
    w = EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'s',R,[0,4,0,0,0,0,0,0,0,0,0,0],ut,0,X0); W{end+1} = sort(w{2},'descend'); fprintf('#');

    % Input Rotations
    w = EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'s',R,[0,0,0,1,0,0,0,0,0,0,0,0],ut,0,X0); W{end+1} = sort(w{2},'descend'); fprintf('#');

    % Extra Input
    w = EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'s',R,[0,0,0,0,0,0,0,1,0,0,0,0],ut,0,X0); W{end+1} = sort(w{2},'descend'); fprintf('#');

    % Parameter Scaling
    w = EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'s',R,[0,0,0,0,0,0,0,0,1,0,0,0],ut,0,X0); W{end+1} = sort(w{2},'descend'); fprintf('#'); % fast drop
    w = EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'s',R,[0,0,0,0,0,0,0,0,2,0,0,0],ut,0,X0); W{end+1} = sort(w{2},'descend'); fprintf('#');

    plotsingvals(W,[2,6,9],'Sensitivity Gramian',14);
    W = {};

    % Baseline
    w = EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'s',R,[0,0,0,0,0,0,0,0,0,1,0,0],ut,0,X0); W{end+1} = sort(w{2},'descend'); fprintf('#');

    % Centering
    w = EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'s',R,[1,0,0,0,0,0,0,0,0,1,0,0],ut,0,X0); W{end+1} = sort(w{2},'descend'); fprintf('#');
    w = EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'s',R,[2,0,0,0,0,0,0,0,0,1,0,0],ut,0,X0); W{end+1} = sort(w{2},'descend'); fprintf('#');
    w = EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'s',R,[3,0,0,0,0,0,0,0,0,1,0,0],ut,0,X0); W{end+1} = sort(w{2},'descend'); fprintf('#');
    w = EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'s',R,[4,0,0,0,0,0,0,0,0,1,0,0],ut,0,X0); W{end+1} = sort(w{2},'descend'); fprintf('#');
    w = EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'s',R,[5,0,0,0,0,0,0,0,0,1,0,0],ut,0,X0); W{end+1} = sort(w{2},'descend'); fprintf('#');
    w = EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'s',R,[6,0,0,0,0,0,0,0,0,1,0,0],ut,0,X0); W{end+1} = sort(w{2},'descend'); fprintf('#');

    % Input scales
    w = EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'s',R,[0,1,0,0,0,0,0,0,0,1,0,0],ut,0,X0); W{end+1} = sort(w{2},'descend'); fprintf('#');
    w = EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'s',R,[0,2,0,0,0,0,0,0,0,1,0,0],ut,0,X0); W{end+1} = sort(w{2},'descend'); fprintf('#');
    w = EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'s',R,[0,3,0,0,0,0,0,0,0,1,0,0],ut,0,X0); W{end+1} = sort(w{2},'descend'); fprintf('#');
    w = EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'s',R,[0,4,0,0,0,0,0,0,0,1,0,0],ut,0,X0); W{end+1} = sort(w{2},'descend'); fprintf('#');

    % State scales
    w = EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'s',R,[0,0,1,0,0,0,0,0,0,1,0,0],ut,0,X0); W{end+1} = sort(w{2},'descend'); fprintf('#');
    w = EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'s',R,[0,0,2,0,0,0,0,0,0,1,0,0],ut,0,X0); W{end+1} = sort(w{2},'descend'); fprintf('#');
    w = EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'s',R,[0,0,3,0,0,0,0,0,0,1,0,0],ut,0,X0); W{end+1} = sort(w{2},'descend'); fprintf('#');
    w = EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'s',R,[0,0,4,0,0,0,0,0,0,1,0,0],ut,0,X0); W{end+1} = sort(w{2},'descend'); fprintf('#');

    % Input Rotations
    w = EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'s',R,[0,0,0,1,0,0,0,0,0,1,0,0],ut,0,X0); W{end+1} = sort(w{2},'descend'); fprintf('#');

    % State Rotations
    w = EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'s',R,[0,0,0,0,1,0,0,0,0,1,0,0],ut,0,X0); W{end+1} = sort(w{2},'descend'); fprintf('#');

    % Extra Input
    w = EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'s',R,[0,0,0,0,0,0,0,1,0,1,0,0],ut,0,X0); W{end+1} = sort(w{2},'descend'); fprintf('#');

    % Parameter Scaling
    w = EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'s',R,[0,0,0,0,0,0,0,0,1,1,0,0],ut,0,X0); W{end+1} = sort(w{2},'descend'); fprintf('#');
    w = EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'s',R,[0,0,0,0,0,0,0,0,2,1,0,0],ut,0,X0); W{end+1} = sort(w{2},'descend'); fprintf('#'); % fast drop

    plotsingvals(W,[2,6,10],'Sensitivity Gramian (Input-Output)',20); fprintf('\n');

%% Identifiability Gramian

    fprintf('Empirical Identifiability Gramian: ');
    W = {};

    % Baseline
    w = EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'i',R,[0,0,0,0,0,0,0,0,0,0,0,0],ut,0,X0); W{end+1} = svd(w{2}); fprintf('#');

    % Centering
    w = EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'i',R,[1,0,0,0,0,0,0,0,0,0,0,0],ut,0,X0); W{end+1} = svd(w{2}); fprintf('#');
    w = EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'i',R,[2,0,0,0,0,0,0,0,0,0,0,0],ut,0,X0); W{end+1} = svd(w{2}); fprintf('#');
    w = EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'i',R,[3,0,0,0,0,0,0,0,0,0,0,0],ut,0,X0); W{end+1} = svd(w{2}); fprintf('#');
    w = EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'i',R,[4,0,0,0,0,0,0,0,0,0,0,0],ut,0,X0); W{end+1} = svd(w{2}); fprintf('#');
    w = EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'i',R,[5,0,0,0,0,0,0,0,0,0,0,0],ut,0,X0); W{end+1} = svd(w{2}); fprintf('#');
    w = EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'i',R,[6,0,0,0,0,0,0,0,0,0,0,0],ut,0,X0); W{end+1} = svd(w{2}); fprintf('#');

    % State scales
    w = EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'i',R,[0,0,1,0,0,0,0,0,0,0,0,0],ut,0,X0); W{end+1} = svd(w{2}); fprintf('#');
    w = EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'i',R,[0,0,2,0,0,0,0,0,0,0,0,0],ut,0,X0); W{end+1} = svd(w{2}); fprintf('#');
    w = EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'i',R,[0,0,3,0,0,0,0,0,0,0,0,0],ut,0,X0); W{end+1} = svd(w{2}); fprintf('#');
    w = EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'i',R,[0,0,4,0,0,0,0,0,0,0,0,0],ut,0,X0); W{end+1} = svd(w{2}); fprintf('#');

    % State Rotations
    w = EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'i',R,[0,0,0,0,1,0,0,0,0,0,0,0],ut,0,X0); W{end+1} = svd(w{2}); fprintf('#');

    % Averaged Identifiability
    w = EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'i',R,[0,0,0,0,0,0,1,0,0,0,0,0],ut,0,X0); W{end+1} = svd(w{2}); fprintf('#');

    % Extra Input
    w = EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'i',R,[0,0,0,0,0,0,0,1,0,0,0,0],ut,0,X0); W{end+1} = svd(w{2}); fprintf('#'); % fast drop

    % Parameter Scaling
    w = EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'i',R,[0,0,0,0,0,0,0,0,1,0,0,0],ut,0,X0); W{end+1} = svd(w{2}); fprintf('#'); % fast drop (r only)
    w = EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'i',R,[0,0,0,0,0,0,0,0,2,0,0,0],ut,0,X0); W{end+1} = svd(w{2}); fprintf('#'); % fast drop (r only)

    % Schur Complement
    w = EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'i',R,[0,0,0,0,0,0,0,0,0,1,0,0],ut,0,X0); W{end+1} = svd(w{2}); fprintf('#');

    plotsingvals(W,[2,6,8],'Identifiability Gramian',13); fprintf('\n');

%% Joint Gramian

    fprintf('Empirical Joint Gramian:           ');
    W = {};

    % Baseline
    w = EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'j',R,[0,0,0,0,0,0,0,0,0,0,0,0],ut,0,X0); W{end+1} = svd(w{2}); fprintf('#');

    % Centering
    w = EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'j',R,[1,0,0,0,0,0,0,0,0,0,0,0],ut,0,X0); W{end+1} = svd(w{2}); fprintf('#');
    w = EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'j',R,[2,0,0,0,0,0,0,0,0,0,0,0],ut,0,X0); W{end+1} = svd(w{2}); fprintf('#');
    w = EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'j',R,[3,0,0,0,0,0,0,0,0,0,0,0],ut,0,X0); W{end+1} = svd(w{2}); fprintf('#');
    w = EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'j',R,[4,0,0,0,0,0,0,0,0,0,0,0],ut,0,X0); W{end+1} = svd(w{2}); fprintf('#');
    w = EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'j',R,[5,0,0,0,0,0,0,0,0,0,0,0],ut,0,X0); W{end+1} = svd(w{2}); fprintf('#');
    w = EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'j',R,[6,0,0,0,0,0,0,0,0,0,0,0],ut,0,X0); W{end+1} = svd(w{2}); fprintf('#');

    % Input scales
    w = EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'j',R,[0,1,0,0,0,0,0,0,0,0,0,0],ut,0,X0); W{end+1} = svd(w{2}); fprintf('#');
    w = EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'j',R,[0,2,0,0,0,0,0,0,0,0,0,0],ut,0,X0); W{end+1} = svd(w{2}); fprintf('#');
    w = EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'j',R,[0,3,0,0,0,0,0,0,0,0,0,0],ut,0,X0); W{end+1} = svd(w{2}); fprintf('#');
    w = EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'j',R,[0,4,0,0,0,0,0,0,0,0,0,0],ut,0,X0); W{end+1} = svd(w{2}); fprintf('#');

    % State scales
    w = EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'j',R,[0,0,1,0,0,0,0,0,0,0,0,0],ut,0,X0); W{end+1} = svd(w{2}); fprintf('#');
    w = EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'j',R,[0,0,2,0,0,0,0,0,0,0,0,0],ut,0,X0); W{end+1} = svd(w{2}); fprintf('#');
    w = EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'j',R,[0,0,3,0,0,0,0,0,0,0,0,0],ut,0,X0); W{end+1} = svd(w{2}); fprintf('#');
    w = EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'j',R,[0,0,4,0,0,0,0,0,0,0,0,0],ut,0,X0); W{end+1} = svd(w{2}); fprintf('#');

    % Input Rotations
    w = EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'j',R,[0,0,0,1,0,0,0,0,0,0,0,0],ut,0,X0); W{end+1} = svd(w{2}); fprintf('#');

    % State Rotations
    w = EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'j',R,[0,0,0,0,1,0,0,0,0,0,0,0],ut,0,X0); W{end+1} = svd(w{2}); fprintf('#');

    % Extra Input
    w = EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'j',R,[0,0,0,0,0,0,0,1,0,0,0,0],ut,0,X0); W{end+1} = svd(w{2}); fprintf('#'); % fast drop

    % Parameter Scaling
    w = EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'j',R,[0,0,0,0,0,0,0,0,1,0,0,0],ut,0,X0); W{end+1} = svd(w{2}); fprintf('#');
    w = EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'j',R,[0,0,0,0,0,0,0,0,2,0,0,0],ut,0,X0); W{end+1} = svd(w{2}); fprintf('#');

    % Schur Complement
    w = EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'j',R,[0,0,0,0,0,0,0,0,0,1,0,0],ut,0,X0); W{end+1} = svd(w{2}); fprintf('#');

    plotsingvals(W,[2,6,11],'Cross-Identifiability Gramian',18);
    W = {};

    % Baseline
    w = EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'j',R,[0,0,0,0,0,0,1,0,0,0,0,0],ut,0,X0); W{end+1} = svd(w{2}); fprintf('#');

    % Centering
    w = EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'j',R,[1,0,0,0,0,0,1,0,0,0,0,0],ut,0,X0); W{end+1} = svd(w{2}); fprintf('#');
    w = EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'j',R,[2,0,0,0,0,0,1,0,0,0,0,0],ut,0,X0); W{end+1} = svd(w{2}); fprintf('#');
    w = EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'j',R,[3,0,0,0,0,0,1,0,0,0,0,0],ut,0,X0); W{end+1} = svd(w{2}); fprintf('#');
    w = EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'j',R,[4,0,0,0,0,0,1,0,0,0,0,0],ut,0,X0); W{end+1} = svd(w{2}); fprintf('#');
    w = EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'j',R,[5,0,0,0,0,0,1,0,0,0,0,0],ut,0,X0); W{end+1} = svd(w{2}); fprintf('#');
    w = EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'j',R,[6,0,0,0,0,0,1,0,0,0,0,0],ut,0,X0); W{end+1} = svd(w{2}); fprintf('#');

    % Input scales
    w = EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'j',R,[0,1,0,0,0,0,1,0,0,0,0,0],ut,0,X0); W{end+1} = svd(w{2}); fprintf('#');
    w = EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'j',R,[0,2,0,0,0,0,1,0,0,0,0,0],ut,0,X0); W{end+1} = svd(w{2}); fprintf('#');
    w = EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'j',R,[0,3,0,0,0,0,1,0,0,0,0,0],ut,0,X0); W{end+1} = svd(w{2}); fprintf('#');
    w = EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'j',R,[0,4,0,0,0,0,1,0,0,0,0,0],ut,0,X0); W{end+1} = svd(w{2}); fprintf('#');

    % State scales
    w = EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'j',R,[0,0,1,0,0,0,1,0,0,0,0,0],ut,0,X0); W{end+1} = svd(w{2}); fprintf('#');
    w = EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'j',R,[0,0,2,0,0,0,1,0,0,0,0,0],ut,0,X0); W{end+1} = svd(w{2}); fprintf('#');
    w = EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'j',R,[0,0,3,0,0,0,1,0,0,0,0,0],ut,0,X0); W{end+1} = svd(w{2}); fprintf('#');
    w = EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'j',R,[0,0,4,0,0,0,1,0,0,0,0,0],ut,0,X0); W{end+1} = svd(w{2}); fprintf('#');

    % Input Rotations
    w = EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'j',R,[0,0,0,1,0,0,1,0,0,0,0,0],ut,0,X0); W{end+1} = svd(w{2}); fprintf('#');

    % State Rotations
    w = EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'j',R,[0,0,0,0,1,0,1,0,0,0,0,0],ut,0,X0); W{end+1} = svd(w{2}); fprintf('#');

    % Extra Input
    w = EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'j',R,[0,0,0,0,0,0,1,1,0,0,0,0],ut,0,X0); W{end+1} = svd(w{2}); fprintf('#'); % fast drop

    % Parameter Scaling
    w = EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'j',R,[0,0,0,0,0,0,1,0,1,0,0,0],ut,0,X0); W{end+1} = svd(w{2}); fprintf('#');
    w = EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'j',R,[0,0,0,0,0,0,1,0,2,0,0,0],ut,0,X0); W{end+1} = svd(w{2}); fprintf('#');

    % Schur Complement
    w = EMGR(LIN,OUT,[M,N,Q],[dt,Tf],'j',R,[0,0,0,0,0,0,1,0,0,1,0,0],ut,0,X0); W{end+1} = svd(w{2}); fprintf('#');

    plotsingvals(W,[2,6,12],'Cross-Identifiability Gramian (Non-symmetric)',18); fprintf('\n\n');
end

function plotsingvals(W,s,l,b)

    if(s(3)==1), figure; end;
    subplot(s(1),s(2),s(3))
    c = winter(numel(W)-1);
    semilogy(W{2}./max(W{2}),'Linewidth',2,'color',c(1,:));
    hold on;
    for k=3:numel(W)

        semilogy(W{k}./max(W{k}),'Linewidth',2,'color',c(k-1,:));
    end
    semilogy(W{1}./max(W{1}),'Linewidth',2,'color',[1,0,0]);
    semilogy(W{b}./max(W{b}),'Linewidth',2,'color',[1,0,1]);
    hold off;
    xlim([1,numel(W{1})]);
    yl = ylim();
    ylim([min(yl(1),0.1),1]);
    xlabel(l)
    set([gca; findall(gca,'Type','text')],'FontSize',6);
end
