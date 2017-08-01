%%% summary: RUNME (Minimal setup and test script)
%%% project: emgr - Empirical Gramian Framework ( gramian.de )
%%% authors: Christian Himpe ( 0000-0003-2194-6754 )
%%% license: 2-Clause BSD (2017)
%$

% Linear System
A = -eye(4)
B = [0;1;0;1]
C = [0,0,1,1]

% (Empircal) Controllability Gramian
emgr(@(x,u,p,t) A*x+B*u,@(x,u,p,t) C*x,[1,4,1],[0.01,1.0],'c')

% (Empirical) Observability Gramian
emgr(@(x,u,p,t) A*x+B*u,@(x,u,p,t) C*x,[1,4,1],[0.01,1.0],'o')

% (Empirical) Cross Gramian
emgr(@(x,u,p,t) A*x+B*u,@(x,u,p,t) C*x,[1,4,1],[0.01,1.0],'x')

% (Empirical) Linear Cross Gramian
emgr(@(x,u,p,t) A*x+B*u,@(x,u,p,t) A'*x+C'*u,[1,4,1],[0.01,1.0],'y')
