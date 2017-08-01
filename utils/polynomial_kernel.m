function w = polynomial_kernel(x,y)
%%% summary: polynomial_kernel - (Polynomial Kernel)
%%% project: emgr - EMpirical GRamian Framework ( http://gramian.de )
%%% authors: Christian Himpe ( 0000-0003-2194-6754 )
%%% license: 2-Clause BSD (2016--2017)
%$
    global EXPONENT;
    if(isempty(EXPONENT)), EXPONENT = 2; end;

    w = (x*y).^EXPONENT + 1.0;
end