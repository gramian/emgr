function w = sigmoid_kernel(x,y)
%%% summary: sigmoid_kernel - (Sigmoid Kernel)
%%% project: emgr - EMpirical GRamian Framework ( http://gramian.de )
%%% authors: Christian Himpe ( 0000-0003-2194-6754 )
%%% license: 2-Clause BSD (2016--2017)
%$
    w = tanh(x*y) + 1.0;
end