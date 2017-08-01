function testall(o)
%%% summary: testall (run all tests and demos)
%%% project: emgr - Empirical Gramian Framework ( gramian.de )
%%% authors: Christian Himpe ( 0000-0003-2194-6754 )
%%% license: 2-Clause BSD (2016--2017)
%$
    if(nargin==0), o = 0; end;

    addpath('tests');
    addpath('demos');
    addpath('utils');

%% Sanity
    disp(' RUNME:');    RUNME;
    disp(' emgrtest:'); emgrtest;

%% Tests
    disp(' test_wx:'); test_wx(o);
    disp(' test_wy:'); test_wy(o);
    disp(' test_wz:'); test_wz(o);
    disp(' test_bt:'); test_bt(o);

    disp(' test_dwx:'); test_dwx(o);
    disp(' test_dwz:'); test_dwz(o);
    disp(' test_wz2:'); test_wz2(o);
    disp(' test_bt2:'); test_bt2(o);

    disp(' test_pwx:'); test_pwx(o);
    disp(' test_pwy:'); test_pwy(o);
    disp(' test_pwz:'); test_pwz(o);
    disp(' test_pbt:'); test_pbt(o);

    disp(' test_ws:');  test_ws(o);
    disp(' test_wi:');  test_wi(o);
    disp(' test_wj:');  test_wj(o);
    disp(' test_wjz:'); test_wjz(o);

    disp(' test_wi2:');  test_wi2(o);
    disp(' test_wj2:');  test_wj2(o);
    disp(' test_cwi:');  test_cwi(o);
    disp(' test_cwj:');  test_cwj(o);

    disp(' test_dwj:');  test_dwj(o);
    disp(' test_kwx:');  test_kwx(o);
    disp(' test_kwy:');  test_kwy(o);
    disp(' test_kpwz:');  test_kpwz(o);

    disp(' test_swx:');  test_swx(o);
    disp(' test_swx2:'); test_swx2(o);
    disp(' test_ewo:');  test_ewo(o);
    disp(' test_ewx:');  test_ewx(o);

%% Demos
    disp(' combined_wj:');   combined_wj(o);
    disp(' gains_wx:');      gains_wx(o);
    disp(' benchmark_ilp:'); benchmark_ilp(o);
    disp(' benchmark_fss:'); benchmark_fss(o);
    disp(' benchmark_lin:'); benchmark_lin(o);
    disp(' benchmark_non:'); benchmark_non(o);
    disp(' energy_wz:');     energy_wz(o);
    disp(' measure:');       measure(o);
    disp(' decentral:');     decentral(o);
    disp(' advection:');     advection(o);
    disp(' nbody:');         nbody(o);
    disp(' blackhole:');     blackhole(o);
end
