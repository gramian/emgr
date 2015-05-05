function y = mysolver(f,g,h,T,x,u,p)
% mysolver (sample custom ode solver for emgr)
% by Christian Himpe, 2014-2015 ( http://gramian.de )
% released under BSD 2-Clause License ( opensource.org/licenses/BSD-2-Clause )
%*

    t = linspace(0,h*T,T);
    U = @(t) u(:,1+min(floor(t/h),T-1));

    % Compute State Trajectory
    if(exist('OCTAVE_VERSION'))
        y = lsode(@(y,t) f(y,U(t),p),x,t);
    else
        [tdummy,y] = ode45(@(t,y) f(y,U(t),p),t,x);
    end;

    y = y';

    if(isnumeric(g) && g==1), return; end;

    % Compute Output Trajectory
    O = numel(g(y(:,1),u(:,1),p));

    for I=1:T
        y(1:O,I) = g(y(:,I),u(:,I),p);
    end;

    y(O+1:end,:) = [];

end

