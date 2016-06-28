function w = trace_kernel(x,y)
% trace_kernel - Kernel for Trace Computation
% Copyright (c) 2016 Christian Himpe ( gramian.de )
% released under BSD 2-Clause License ( opensource.org/licenses/BSD-2-Clause )
%*
    w = sum(sum(x.*y'));
end