function blackhole(o)
% stable orbit parameter identification inside event horizon of black hole
% by Christian Himpe, 2013 ( http://gramian.de )
% released under BSD 2-Clause License ( opensource.org/licenses/BSD-2-Clause )
%*

if(exist('emgr')~=2) disp('emgr framework is required. Download at http://gramian.de/emgr.m'); return; end

%%%%%%%% Setup %%%%%%%%

 t = [0 0.005 5];
 T = (t(3)-t(1))/t(2);
 u = zeros(1,T);
 R = 14;

% PLANET
 xu = [0.4;pi/2;0;0];
 pu = [0.568;1.13;0.13;0.9982;0.05;0;1]; %E,L,Q,a,e,eps,mu

% PHOTON
 xp = [0.2;pi/2;0;0];
 EE = 0.568; %10.5
 pp = [EE;1.38*EE;0.03*EE*EE;0.9982;0.05;0;0]; %E,L,Q,a,e,eps,mu

%%%%%%%% Identification %%%%%%%%

% FULL
 Y = [rk2(@orbit,@bl2c,[0 4 3],t,xu,u,pu);rk2(@orbit,@bl2c,[0 4 3],t,xp,u,pp)];

% PLANET
 WS = emgr(@orbit,@bl2c,[0 4 3],pu,t,'s',[1 0 0 0 0 0 0 0 0 0],1,0,xu);
 PLANET_SENSITIVITY = diag(WS{2})

% PHOTON
 WS = emgr(@orbit,@bl2c,[0 4 3],pp,t,'s',[1 0 0 0 0 0 0 0 0 0],1,0,xp);
 PHOTON_SENSITIVITY = diag(WS{2})

%%%%%%%% Output %%%%%%%%

 if(nargin<1 || o==0 ) return; end

 figure('PaperSize',[7,7],'PaperPosition',[0,0,7,7]);
 grid on;
 hold on;
 p0 = plot3(0,0,0,'*','Color','black');				%singularity
 p1 = plot3(Y(1,end),Y(2,end),Y(3,end),'*','Color','red');		%planet
 p2 = plot3(Y(1,:),Y(2,:),Y(3,:),'Color','red');			%planet orbit
 p3 = plot3(Y(4,end),Y(5,end),Y(6,end),'*','Color','blue');		%photon
 p4 = plot3(Y(4,:),Y(5,:),Y(6,:),'Color','blue');			%photon orbit
 lg = legend([p0 p2 p4],'singularity','planet orbit','photon orbit'); set(lg,'FontSize',5);
 hold off;
 xl = ceil(10*max([abs(Y(1,:)),abs(Y(4,:))]))*0.1;
 yl = ceil(10*max([abs(Y(2,:)),abs(Y(5,:))]))*0.1;
 zl = ceil(10*max([abs(Y(3,:)),abs(Y(6,:))]))*0.1;
 set(gca,'Xlim',[-xl xl],'Ylim',[-yl yl],'Zlim',[-zl zl]);

 if(o==2 && ~exist('OCTAVE_VERSION')), camorbit(-30,-60); print -dpng blackhole.png; end

%%%%%%%% Orbit %%%%%%%%

function x = orbit(x,u,p)

 E = p(1);
 L = p(2);
 Q = p(3);
 a = p(4);
 e = p(5);
 w = p(6);
 m = p(7);

 D  = x(1)^2 - 2*x(1) + a^2 + e^2;
 S  = x(1)^2 + a^2*cos(x(2))^2;
 P  = E*(x(1)^2 + a^2) + e*w*x(1) - a*L;
 Vt = Q - cos(x(2))^2*(a^2*(m^2 - E^2) + L^2*sin(x(2))^(-2) );
 Vr = P^2 - D*(m^2*x(1)^2 + (L - a*E)^2 + Q);

 x = abs([ sqrt(Vr) ; sqrt(Vt) ; L*sin(x(2))^(-2)+a*(P/D-E) ; a*(L-a*E*sin(x(2))^2)+P/D*(x(1)^2+a^2) ]./S);

%%%%%%%% Boyer-Lindquist to Cartesian %%%%%%%%

function y = bl2c(x,u,p)

 a = p(4);

 y = [ sqrt(x(1)^2+a^2)*sin(x(2))*cos(x(3)) ; sqrt(x(1)^2+a^2)*sin(x(2))*sin(x(3)) ; x(1)*cos(x(2)) ];

%%%%%%%% Integrator %%%%%%%%

function y = rk2(f,g,q,t,x,u,p)

 T = (t(3)-t(1))/t(2);
 y = zeros(q(3),T);
 h = t(2);

 for A=1:T
  x = x + h*f(x + 0.5*h*f(x,u(:,A),p),u(:,A),p); %Improved Eulers Method
  y(:,A) = g(x,u(:,A),p);
 end
