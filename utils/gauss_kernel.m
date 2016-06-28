function w = gauss_kernel(x,y)
% gaussian_kernel - Gaussian Kernel
% Copyright (c) 2016 Christian Himpe ( gramian.de )
% released under BSD 2-Clause License ( opensource.org/licenses/BSD-2-Clause )
%*
    global GAMMA;
    if(isempty(GAMMA)), GAMMA = 1; end;

    r = sqrt(1.0/GAMMA)*(x - y');
    w = exp(-(r*r').^2.0);
end