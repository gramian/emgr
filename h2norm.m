function e = h2norm(A,B,C,h,T)
% h2norm (H2 norm of linear control system)
% by Christian Himpe, 2014-2015 ( http://gramian.de )
% released under BSD 2-Clause License ( opensource.org/licenses/BSD-2-Clause )
%*
    if(exist('emgr')~=2)
        disp('emgr framework is required. Download at http://gramian.de/emgr.m');
        return;
    end

    WC = emgr(@(x,u,p) A*x+B*u,1,[size(B,2) size(A,1) size(C,1)],[0 h T],'c');
    CC = C'*C;

    e = sqrt(sum(sum(CC.*WC)));
end
