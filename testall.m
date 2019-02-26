%%% summary: testall
%%% project: emgr - Empirical Gramian Framework ( https://gramian.de )
%%% authors: Christian Himpe ( 0000-0003-2194-6754 )
%%% license: BSD-2-Clause (2019)

if(exist('OCTAVE_VERSION','builtin'))
    svd_driver('gesvd');
end

RUNME

emgrtest();

moretests('i');
moretests('s');
moretests('c');
moretests('r');

examples('hnm');
examples('isp');
examples('fss');
examples('nrc');
examples('rqo');
examples('lte');
examples('fbc');
examples('qso');
