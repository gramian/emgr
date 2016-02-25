function blackhole(o)
% stable orbit parameter identification inside the event horizon of black hole
% by Christian Himpe, 2013-2016 ( http://gramian.de )
% released under BSD 2-Clause License ( opensource.org/licenses/BSD-2-Clause )
%*
    if(exist('emgr')~=2)
        error('emgr not found! Get emgr at: http://gramian.de');
    else
        global ODE; ODE = [];
        fprintf('emgr (version: %g)\n',emgr('version'));
    end

%% SETUP
    T = [0.005,5.0];
    L = floor(T(2)/T(1)) + 1;
    U = zeros(1,L);

   % Planet Parameters
    xu = [0.4;pi/2;0;0];
    pu = [0.568;1.13;0.13;0.9982;0.05;0;1]; % E,L,Q,a,e,eps,m

    % Photon Parameters
    xp = [0.2;pi/2;0;0];
    EE = 0.568; %10.5;
    pp = [EE;1.38*EE;0.03*EE*EE;0.9982;0.05;0;0]; % E,L,Q,a,e,eps,m

%% MAIN
    Y = [ODE(@orbit,@bl2c,T,xu,U,pu);... % Full Order Planet
         ODE(@orbit,@bl2c,T,xp,U,pp)];   % Full Order Photon

    fprintf('Parameters: E,L,Q,a,e,eps,m\n');

    % Planet Sensitivities
    WS = emgr(@orbit,@bl2c,[0,4,3],T,'s',[0.9*pu,1.1*pu],[1,0,0,0,0,0,0,0,0,0,0,0],1,0,xu);
    PLANET_SENSITIVITY = full(diag(WS{2}))

    % Photon Sensitivities
    WS = emgr(@orbit,@bl2c,[0,4,3],T,'s',[0.9*pp,1.1*pp],[1,0,0,0,0,0,0,0,0,0,0,0],1,0,xp);
    PHOTON_SENSITIVITY = full(diag(WS{2}))

    %% Output
    if(nargin>0 && o==0), return; end; 
    figure('Name',mfilename,'NumberTitle','off');
    grid on;
    hold on;
    p0 = plot3(0,0,0,'*','Color','black');                     %singularity
    p1 = plot3(Y(1,end),Y(2,end),Y(3,end),'*','Color','red');  %planet
    p2 = plot3(Y(1,:),Y(2,:),Y(3,:),'Color','red');            %planet orbit
    p3 = plot3(Y(4,end),Y(5,end),Y(6,end),'*','Color','blue'); %photon
    p4 = plot3(Y(4,:),Y(5,:),Y(6,:),'Color','blue');           %photon orbit
    l = legend([p0 p2 p4],'singularity','planet orbit','photon orbit');
    set(l,'FontSize',5);
    hold off;
    xl = ceil(10*max([abs(Y(1,:)),abs(Y(4,:))]))*0.1;
    yl = ceil(10*max([abs(Y(2,:)),abs(Y(5,:))]))*0.1;
    zl = ceil(10*max([abs(Y(3,:)),abs(Y(6,:))]))*0.1;
    set(gca,'Xlim',[-xl,xl],'Ylim',[-yl,yl],'Zlim',[-zl,zl]);
    view(-30,30);
    if(nargin>0 && o==1), print('-dsvg',[mfilename(),'.svg']); end;
end

%% ======== Orbit ========
function x = orbit(x,u,p)

    E = p(1); % E
    L = p(2); % L
    Q = p(3); % Q
    a = p(4); % a
    e = p(5); % e
    w = p(6); % m
    m = p(7); % eps

    D  = x(1)^2 - 2*x(1) + a^2 + e^2;
    S  = x(1)^2 + a^2*cos(x(2))^2;
    P  = E*(x(1)^2 + a^2) + e*w*x(1) - a*L;
    Vt = Q - cos(x(2))^2*(a^2*(m^2 - E^2) + L^2*sin(x(2))^(-2) );
    Vr = P^2 - D*(m^2*x(1)^2 + (L - a*E)^2 + Q);

    x = abs([ sqrt(Vr) ; ... 
              sqrt(Vt) ; ... 
              L*sin(x(2))^(-2)+a*(P/D-E) ; ...
              a*(L-a*E*sin(x(2))^2)+P/D*(x(1)^2+a^2) ]./S);
end

%% ======== Boyer-Lindquist to Cartesian ========
function y = bl2c(x,u,p)

    a = p(4);

    y = [ sqrt(x(1)^2+a^2)*sin(x(2))*cos(x(3)) ; ...
          sqrt(x(1)^2+a^2)*sin(x(2))*sin(x(3)) ; ...
          x(1)*cos(x(2)) ];
end
