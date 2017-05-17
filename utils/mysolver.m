function y = mysolver(f,g,t,x0,u,p)
%%% summary: mysolver (sample custom solver for emgr)
%%% project: emgr - Empirical Gramian Framework ( http://gramian.de )
%%% authors: Christian Himpe ( 0000-0003-2194-6754 )
%%% license: 2-Clause BSD (2016--2017)
%$
    h = t(1);
    T = t(2);
    L = floor(T/h) + 1;

    % Compute State Trajectory
    [k,x] = ode45(@(t,x) f(x,u(t),p,t),[0,T],x0)
    x = x';

    if(isnumeric(g) && g==1),
        z = x;
    else % Compute Output Trajectory
        z = g(x(:,1),u(:,1),p,0);
        K = numel(k)
        z(end,K) = 0;
        for I=2:K
            z(:,I) = g(x(:,I),u(t),p,k);
        end;
    end
    y = interp1(k,z',0:h:T)';
end
