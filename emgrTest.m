%%% project: emgr - EMpirical GRamian Framework ( https://gramian.de )
%%% version: 5.99 (2022-04-13)
%%% authors: Christian Himpe (0000-0003-2194-6754)
%%% license: BSD-2-Clause (opensource.org/licenses/BSD-2-Clause)
%%% summary: emgrTest - run all system tests

tid = tic();

RUNME;

disp(' ');

estTest('matrix_equation');
estTest('singular_values');
estTest('model_reduction');
estTest('parameter_reduction');
estTest('combined_reduction');
estTest('decentralized_control');
estTest('state_sensitivity');
estTest('parameter_sensitivity');
estTest('parameter_identifiability');
estTest('uncertainty_quantification');
estTest('nonlinearity_quantification');
estTest('gramian_index');
estTest('system_index');
estTest('system_norm');
estTest('tau_function');

disp(' ');

estDemo('hnm');
estDemo('isp');
estDemo('fss');
estDemo('nrc');
estDemo('rqo');
estDemo('aps');
estDemo('lte');
estDemo('fbc');
estDemo('qso');

disp(' ');

%estProbe([],'controllability','linear')
%estProbe([],'minimality','linear')
%estProbe([],'observability','nonlinear')
%estProbe([],'minimality','nonlinear')

totalTime = toc(tid)

