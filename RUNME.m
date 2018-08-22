%%% summary: RUNME (Minimal test script)
%%% project: emgr - Empirical Gramian Framework ( https://gramian.de )
%%% authors: Christian Himpe ( 0000-0003-2194-6754 )
%%% license: 2-Clause BSD (2017--2018)
%$

if(exist('OCTAVE_VERSION','builtin'))
    EMGR = @emgr_oct;
elseif(verLessThan('matlab','9.1'))
    EMGR = @emgr_lgc;
else
    EMGR = @emgr;
end

emgr_version = emgr('version')

% Linear System
A = -eye(4)
B = [0;1;0;1]
C = [0,0,1,1]

% (Empircal) Controllability Gramian
WC = EMGR(@(x,u,p,t) A*x+B*u,@(x,u,p,t) C*x,[1,4,1],[0.01,1.0],'c')

% (Empirical) Observability Gramian
WO = EMGR(@(x,u,p,t) A*x+B*u,@(x,u,p,t) C*x,[1,4,1],[0.01,1.0],'o')

% (Empirical) Cross Gramian
WX = EMGR(@(x,u,p,t) A*x+B*u,@(x,u,p,t) C*x,[1,4,1],[0.01,1.0],'x')

% (Empirical) Linear Cross Gramian
WY = EMGR(@(x,u,p,t) A*x+B*u,@(x,u,p,t) A'*x+C'*u,[1,4,1],[0.01,1.0],'y')
