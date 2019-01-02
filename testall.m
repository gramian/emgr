%%% summary: testall
%%% project: emgr - Empirical Gramian Framework ( https://gramian.de )
%%% authors: Christian Himpe ( 0000-0003-2194-6754 )
%%% license: BSD-2-Clause (2019)

emgrtest();

moretests(1);
moretests(0);
moretests(Inf);

examples('hnm');
examples('isp');
examples('fss');
examples('nrc');
examples('lte');
examples('fbc');
examples('qso');
