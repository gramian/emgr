%%% project: emgr - EMpirical GRamian Framework ( https://gramian.de )
%%% version: 5.99 (2022-04-13)
%%% authors: Christian Himpe (0000-0003-2194-6754)
%%% license: BSD-2-Clause (opensource.org/licenses/BSD-2-Clause)
%%% summary: RUNME - minimal test script

emgr_VERSION = emgr('version')

% Linear System
A = -eye(4)					% System Matrix
B = [0;1;0;1]					% Input Matrix
C = [0,0,1,1]					% Output Matrix

P = zeros(4,1);				% Parameter Vector
Q = [0.01*ones(4,1),[0.25;0.5;0.75;1.0]];	% Min and Max Parameter Box

f = @(x,u,p,t) A*x + B*u + p;			% (Affine) Linear Vector Field
g = @(x,u,p,t) C*x;				% Linear Output Functional
h = @(x,u,p,t) A'*x + C'*u;			% Adjoint Vector Field

s = [1,4,1];					% System dimension
t = [0.01,1.0];				% Time discretization

% (Empircal) Controllability Gramian
WC = emgr(f,g,s,t,'c',P)

% (Empirical) Observability Gramian
WO = emgr(f,g,s,t,'o',P)

% (Empirical) Cross Gramian
WX = emgr(f,g,s,t,'x',P)

% (Empirical) Linear Cross Gramian
WY = emgr(f,h,s,t,'y',P)

% (Empirical) Sensitivity Gramian
WCWS = emgr(f,g,s,t,'s',Q)

% (Empirical) Identifiability Gramian
WOWI = emgr(f,g,s,t,'i',Q)

% (Empirical) Cross-Identifiability Gramian
WXWJ = emgr(f,g,s,t,'j',Q)

