function fulltest(o)
%%% summary: fulltest (run all tests and demos)
%%% project: emgr - Empirical Gramian Framework ( gramian.de )
%%% authors: Christian Himpe ( 0000-0003-2194-6754 )
%%% license: 2-Clause BSD (2016)
%$
    if(nargin==0), o = 0; end;

    addpath('tests');
    addpath('demos');
    addpath('utils');

%% Sanity
    emgrtest;

%% Tests
    test_wx(o);
    test_wy(o);
    test_bt(o);
    test_bt2(o);

    test_wz(o);
    test_wz2(o);
    test_kwx(o);
    test_swx(o)

    test_pwx(o);
    test_pwy(o);
    test_pbt(o);
    test_pwz(o);

    test_ws(o);
    test_wi(o);
    test_wj(o);
    test_wjz(o)

    test_ws2(o);
    test_wi2(o);
    test_wj2(o);
    test_wjz2(o);

    test_cws(o);
    test_cwi(o);
    test_cwj(o);
    test_cwjz(o);

    test_dwx(o);
    test_dwz(o);
    test_dwj(o);
    test_dwjz(o);

%% Demos

    combined_wj(o);
    gains_wx(o);
    benchmark_ilp(o);
    benchmark_lin(o);
    benchmark_non(o);
    energy_wz(o);
    measure(o);
    decentral(o);
    advection(o);
    nbody(o);
    blackhole(o);
end
