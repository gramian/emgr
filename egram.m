function W = egram(a,b,c,w,varargin)
% egram (shortcut encapsulation)
% by Christian Himpe, 2013-2014 ( http://gramian.de )
% released under BSD 2-Clause License ( opensource.org/licenses/BSD-2-Clause )
%*

if(exist('emgr')~=2) disp('emgr framework is required. Download at http://gramian.de/emgr.m'); return; end

J = size(b,2);
N = size(a,1);
O = size(c,1);

if( (N~=size(a,2)) || (N~=size(b,1)) || (N~=size(c,2)) ), error('ERROR: matrix dimension mismatch'); end;

W = emgr(@(x,u,p) a*x+b*u,@(x,u,p) c*x,[J N O],[0 (1.0/N) 1.0],w,varargin);
