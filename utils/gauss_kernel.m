function w = gauss_kernel(x,y)
%%% summary: gauss_kernel - (Gaussian Kernel)
%%% project: emgr - EMpirical GRamian Framework ( http://gramian.de )
%%% authors: Christian Himpe ( 0000-0003-2194-6754 )
%%% license: 2-Clause BSD (2016--2017)
%$
    global GAMMA;
    if(isempty(GAMMA)), GAMMA = 1; end;

    r = sqrt(1.0/GAMMA)*(x - y');
    w = exp(-(r*r'));
end