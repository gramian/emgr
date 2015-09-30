function vernval(J)
% vernval (verification and validation)
% by Christian Himpe, 2013-2015 ( http://gramian.de )
% released under BSD 2-Clause License ( opensource.org/licenses/BSD-2-Clause )
%*
    if(exist('emgr')~=2)
        disp('emgr framework is required. Download at http://gramian.de/emgr.m');
        return;
    end

    rand('seed',1009);

    ok = @(t) fprintf([t,' - OK \n']);

    if(nargin<1) J=4; end	% number of inputs
    N = J*J;			% number of states
    O = J;			% number of outputs

    S = 0;			% start time
    h = 0.01;			% time step
    T = 1.0;			% end time

    A = rand(N,N);		% random system matrix
    A(1:N+1:end) = -0.55*N;	% ensure stability
    A = 0.5*(A+A');		% symmetrize system matrix
    B = rand(N,J);		% random input matrix
    C = B';			% ensure state-space symmetric system
    P = ones(N,1);		% parameter vector

    f = @(x,u,p) A*x+B*u;   	% linear dynamic system vector field
    g = @(x,u,p) C*x;       	% linear output functional

    G = @(x,u,p) A'*x+C'*u;	% adjoint dynamic system vector field

    F = @(x,u,p) A*x+B*u+p; 	% linear parametrized system vector field


    WC = emgr(f,g,[J,N,O],[S,h,T],'c',0,1); ok('Empirical Controllability Gramian: WC');
    Wc = emgr(G,g,[J,N,O],[S,h,T],'c',0,1); ok('Adjoint Empirical Controllability Gramian: Wc');
    WO = emgr(f,g,[J,N,O],[S,h,T],'o',0,1); ok('Empirical Observability Gramian: WO');
    WX = emgr(f,g,[J,N,O],[S,h,T],'x',0,1); ok('Empirical Cross Gramian: WX');
    WY = emgr(f,G,[J,N,O],[S,h,T],'y',0,1); ok('Empirical Linear Cross Gramian: WY');

    dWCWc = norm(WC-Wc,'fro')
    dWCWY = norm(WC-WY,'fro')
    dWCWO = norm(WC-WO,'fro')
    dWCWX = norm(WC-WX,'fro')
    dWOWX = norm(WO-WX,'fro')
    dWXWY = norm(WX-WY,'fro')

    WC = emgr(F,g,[J,N,O],[S,h,T],'c',P); ok('Parametrized Empirical Controllability Gramian: WC');
    WO = emgr(F,g,[J,N,O],[S,h,T],'o',P); ok('Parametrized Empirical Observability Gramian: WO');
    WX = emgr(F,g,[J,N,O],[S,h,T],'x',P); ok('Parametrized Empirical Cross Gramian: WX');
    WS = emgr(F,g,[J,N,O],[S,h,T],'s',P); ok('Empirical Sensitivity Gramian: WS');
    WI = emgr(F,g,[J,N,O],[S,h,T],'i',P); ok('Empirical Identifiability Gramian: WI');
    WJ = emgr(F,g,[J,N,O],[S,h,T],'j',P); ok('Empirical Joint Gramian: WJ - OK');

    dWCWS = norm(WC-WS{1},'fro')
    dWOWI = norm(WO-WI{1},'fro')
    dWXWJ = norm(WX-WJ{1},'fro')
    dWIWJ = norm(WI{2}-WJ{2},'fro')
end
