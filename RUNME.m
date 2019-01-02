%%% summary: RUNME (Minimal test script)
%%% project: emgr - Empirical Gramian Framework ( https://gramian.de )
%%% authors: Christian Himpe ( 0000-0003-2194-6754 )
%%% license: BSD-2-Clause (2019)

if(exist('OCTAVE_VERSION','builtin'))	% Use Octave Version
    EMGR = @emgr_oct;
elseif(verLessThan('matlab','9.1'))	% Use Matlab Legacy Version
    EMGR = @emgr_lgc;
else					% Use Standard Version
    EMGR = @emgr;
end

EMGR_VERSION = EMGR('version')

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
WC = EMGR(f,g,[1,4,1],[0.01,1.0],'c',P)

% (Empirical) Observability Gramian
WO = EMGR(f,g,[1,4,1],[0.01,1.0],'o',P)

% (Empirical) Cross Gramian
WX = EMGR(f,g,[1,4,1],[0.01,1.0],'x',P)

% (Empirical) Linear Cross Gramian
WY = EMGR(f,h,[1,4,1],[0.01,1.0],'y',P)

% (Empirical) Sensitivity Gramian
WCWS = EMGR(f,g,[1,4,1],[0.01,1.0],'s',Q); WS = WCWS{2}

% (Empirical) Identifiability Gramian
WOWI = EMGR(f,g,[1,4,1],[0.01,1.0],'i',Q); WI = WOWI{2}

% (Empirical) Cross-Identifiability Gramian
WXWJ = EMGR(f,g,[1,4,1],[0.01,1.0],'j',Q); WJ = WXWJ{2}
