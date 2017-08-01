function y = mysolver(f,g,t,x0,u,p)
%%% summary: mysolver (sample custom solver for emgr)
%%% project: emgr - EMpirical GRamian Framework ( http://gramian.de )
%%% authors: Christian Himpe ( 0000-0003-2194-6754 )
%%% license: 2-Clause BSD (2016--2017)
%$
    h = t(1);
    T = t(2);

    % Compute State Trajectory
    [S,x] = ode45(@(t,x) f(x,u(t),p,t),[h,T],x0);

    % Compute Output Trajectory
    K = numel(S);
    z = g(x(1,:)',u(S(1)),p,S(1));
    z(end,K) = 0;
    for k=2:K
        tk = S(k);
        z(:,k) = g(x(k,:)',u(tk),p,tk);
    end;

    y = interp1(S,z',0:h:T)';
end
