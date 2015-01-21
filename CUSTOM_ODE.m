function y = CUSTOM_ODE(f,h,T,x,u,p)
% CUSTOM_ODE (sample custom ode solver for emgr)
% by Christian Himpe, 2014-2015 ( http://gramian.de )
% released under BSD 2-Clause License ( opensource.org/licenses/BSD-2-Clause )
%*

 H = 1.0./h;
 t = linspace(0,h*T,T);

 if(exist('OCTAVE_VERSION'))
  y = lsode(@(y,t) f(y,u(:,1+min(round(t*H),T-1)),p),x,t)';
 else,
  y = deval(ode45(@(t,y) f(y,u(:,1+min(round(t*H),T-1)),p),[0 h*T],x),t);
 end;
