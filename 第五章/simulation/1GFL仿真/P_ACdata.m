function [baseMVA, bus, gen, branch] = ACdata
%case 5 nodes    Power flow data for 5 bus, 2 generator case.
%   Please see 'help caseformat' for details on the case file format.
%
%   case file can be used together with dc case files "case5_stagg_....m"
%
%   Network data from ...
%   G.W. Stagg, A.H. El-Abiad, "Computer methods in power system analysis",
%   McGraw-Hill, 1968.
%
%   MATPOWER case file data provided by Jef Beerten.

%% MATPOWER Case Format : Version 1
%%-----  Power Flow Data  -----%%
%% system MVA base
baseMVA = 1000;

%% bus data
%  1        2       3       4       5   6       7       8       9   10         11       12  13
%  BUS_I, BUS_TYPE, PD,     QD,     GS, BS,  BUS_AREA,  VM,     VA, BASE_KV, ZONE,    VMAX, VMIN
%	bus_i	type	Pd      Qd      Gs	Bs      area	Vm      Va	baseKV	 zone     Vmax	Vmin
bus = [
	1       3       0       0       0	0       1       1.0     0      0       1       1.1     0.9;
	2       1       0       0   	0	0       1       1.0 	 0      0       1       1.1     0.9;
% 	3       1       0   	0       0	0       1       1.0 	 0   	0       1       1.1     0.9;
% 	4       2       0   	0   	0	0       1       1.05	 0      0       1       1.1     0.9;
% 	5       1       0       0       0	0       1       1.00	 0    	0       1       1.1     0.9;
% 	6       1       0       0    	0	0       1       1.00	 0      0       1       1.1     0.9;
% 	7       1       967     -100    0	0       1       1.00	 0      0       1       1.1 	0.9;
% 	8       1       0       0       0	0       1       1.00	 0      0       1       1.1 	0.9;   
%     9       1       1767  	-250 	0	0       1       1.00	 0      0       1       1.1     0.9;
% 	10      1       0       0       0	0       1       1.00	 0      0       1       1.1 	0.9;
% 	11      1       0       0       0	0       1       1.00	 0      0       1       1.1     0.9;	
  
    ];
%-400MVar£»250MW£»120MW£¬140Mvar£»
%% generator data
%   1       2       3       4      5    6       7           8           9   10
% GEN_BUS, PG,    QG,      QMAX, QMIN, VG,      MBASE,    GEN_STATUS, PMAX, PMIN,
%	bus	  Pg      Qg	   Qmax	 Qmin	Vg      mBase       status	Pmax	Pmin
gen = [
    1    0     0   1000     -400    1.0     100     1       1000     0;
    2    1000     0    5000     -400    1.0     100     1       1000 	 0;
%     3    0     0    4000     -400    1.03     100 	1       1000     0;
%     4    650     0    2400     -400    1.01     100 	1       1000     0;
    ];

%% branch data
%  1        2       3       4           5           6       7       8       9       10      11          12       13
%  F_BUS,  T_BUS, BR_R,    BR_X,      BR_B,        RATE_A, RATE_B,  RATE_C, TAP,   SHIFT,  BR_STATUS, ANGMIN, ANGMAX,
%	fbus	tbus	r       x           b           rateA	rateB	rateC	ratio	angle	status
branch = [
    1       2       0.0   0.3454        0.1      0       0         0      0      0      1;
%     2       3       0.0   0.25        0.1      0       0         0      0      0      1;
%     3       11      0       0.01670     0           0       0         0      1   	0       1;
%     4       10      0       0.01670     0           0       0         0      1   	0       1;
%     5       6       0.0025	0.025       0.04375     0       0         0      0   	0       1;
%     6       7       0.001	0.01        0.0175      0       0         0      0   	0       1;
%     7       8       0.011   0.11        0.1925      0       0         0      0      0      1;
%     7       8       0.011   0.11        0.1925      0       0         0      0      0      1;
%     8       9       0.011   0.11        0.1925      0       0         0      0      0      1;
%     8       9       0.011   0.11        0.1925      0       0         0      0      0      1;
%     9       10      0.001	0.01        0.0175      0       0         0      0   	0       1;
%     10      11      0.0025	0.025       0.04375     0       0         0      0   	0       1;
    1       1       inf   inf        0.1      0       0         0      0      0      1;
    ];
return;
