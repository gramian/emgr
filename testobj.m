function testobj()
%%% summary: testobj (test state-space obj)
%%% project: emgr - Empirical Gramian Framework ( https://gramian.de )
%%% authors: Christian Himpe ( 0000-0003-2194-6754 )
%%% license: BSD-2-Clause (2019)  

    N = 16;

    A = spdiags(-rand(N,1),0,N,N);
    B = ones(N,1);
    C = B';
    D = 0;
    sys = ss(A,B,C,D);

    curios(sys,'state-reduction','linear-balanced-truncation',{'noscore'});

    sys = ss(A,B,C,D);
    sys.e = speye(N);

    curios(sys,'state-reduction','linear-balanced-truncation',{'noscore'});
end