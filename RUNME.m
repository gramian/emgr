%%% project: emgr - EMpirical GRamian Framework ( https://gramian.de )
%%% version: 5.8 (2020-05-01)
%%% authors: Christian Himpe (0000-0003-2194-6754)
%%% license: BSD-2-Clause (opensource.org/licenses/BSD-2-Clause)
%%% summary: RUNME - minimal test script

emgr_VERSION = emgr('version')

% Linear System
A = -eye(4)				% System Matrix
B = [0;1;0;1]				% Input Matrix
C = [0,0,1,1]				% Output Matrix

P = zeros(4,1);				% Parameter Vector
Q = [zeros(4,1),[0.25;0.5;0.75;1.0]];	% Min and Max Parameter Box

f = @(x,u,p,t) A*x + B*u + p;		% (Affine) Linear Vector Field
g = @(x,u,p,t) C*x;			% Linear Output Functional
h = @(x,u,p,t) A'*x + C'*u;		% Adjoint Vector Field

% (Empircal) Controllability Gramian
WC = emgr(f,g,[1,4,1],[0.01,1.0],'c',P)

% (Empirical) Observability Gramian
WO = emgr(f,g,[1,4,1],[0.01,1.0],'o',P)

% (Empirical) Cross Gramian
WX = emgr(f,g,[1,4,1],[0.01,1.0],'x',P)

% (Empirical) Linear Cross Gramian
WY = emgr(f,h,[1,4,1],[0.01,1.0],'y',P)

% (Empirical) Sensitivity Gramian
WCWS = emgr(f,g,[1,4,1],[0.01,1.0],'s',Q); WS = WCWS{2}

% (Empirical) Identifiability Gramian
WOWI = emgr(f,g,[1,4,1],[0.01,1.0],'i',Q); WI = WOWI{2}

% (Empirical) Cross-Identifiability Gramian
WXWJ = emgr(f,g,[1,4,1],[0.01,1.0],'j',Q); WJ = WXWJ{2}
