function w = trace_kernel(x,y)
%%% summary: trace_kernel (Pseudo-Kernel for Trace Computation)
%%% project: emgr - Empirical Gramian Framework ( http://gramian.de )
%%% authors: Christian Himpe ( 0000-0003-2194-6754 )
%%% license: 2-Clause BSD (2016)
%$
    w = sum(sum(x.*y'));
end