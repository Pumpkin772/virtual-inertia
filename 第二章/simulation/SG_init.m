function [GSGInit] = SG_init(resultsac,Syn_Machidx_g)


%%%  Bus No. Xd  Xd1  Xq  Td1 J D
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SGdata=...
[
    30   1.81*0.1  0.3*0.1   1.78*0.1 8 4.5 5;
    31   1.81*0.1  0.3*0.1   1.78*0.1 8 11.6 5;
    32   1.81*0.1  0.3*0.1   1.78*0.1 8 11.6 5;
    33   1.81*0.1  0.3*0.1   1.78*0.1 8 4.5 5;
    34   1.81*0.1  0.3*0.1   1.78*0.1 8 11.6 5;
    35   1.81*0.1  0.3*0.1   1.78*0.1 8 1.9 5;
    36   1.81*0.1  0.3*0.1   1.78*0.1 8 1.7 5;
    37   1.81*0.1  0.3*0.1   1.78*0.1 8 11.6 5;
    38   1.81*0.1  0.3*0.1   1.78*0.1 8 11.6 5;
    39   1.81*0.1  0.3*0.1   1.78*0.1 8 11.6 5;
];

Xd = SGdata(:,2);
Xd1 = SGdata(:,3);
Xq = SGdata(:,4);
Td1 = SGdata(:,5);
J  = SGdata(:,6);
D  = SGdata(:,7);

GSGInit.Xd = Xd;
GSGInit.Xd1 = Xd1;
GSGInit.Xq = Xq;
GSGInit.Td1 = Td1;
GSGInit.J  = J;
GSGInit.D  = D;

%% 
%  Ksg    K1    K3    K5    K7    T1   T2    T3    T4  T5    T6    T7 KAVR TAVR
governordata=...
[
20.00  0.22  0.22  0.30  0.26  0.2  0.00  0.10  0.25  4.00  4.00  0.40 200 0.015;
20.00  0.22  0.22  0.30  0.26  0.2  0.00  0.10  0.25  4.00  4.00  0.40 200 0.015;
20.00  0.22  0.22  0.30  0.26  0.2  0.00  0.10  0.25  4.00  4.00  0.40 200 0.015;
20.00  0.22  0.22  0.30  0.26  0.2  0.00  0.10  0.25  4.00  4.00  0.40 200 0.015;
20.00  0.22  0.22  0.30  0.26  0.2  0.00  0.10  0.25  4.00  4.00  0.40 200 0.015;
20.00  0.22  0.22  0.30  0.26  0.2  0.00  0.10  0.25  4.00  4.00  0.40 200 0.015;
20.00  0.22  0.22  0.30  0.26  0.2  0.00  0.10  0.25  4.00  4.00  0.40 200 0.015;
20.00  0.22  0.22  0.30  0.26  0.2  0.00  0.10  0.25  4.00  4.00  0.40 200 0.015;
20.00  0.22  0.22  0.30  0.26  0.2  0.00  0.10  0.25  4.00  4.00  0.40 200 0.015;
20.00  0.22  0.22  0.30  0.26  0.2  0.00  0.10  0.25  4.00  4.00  0.40 200 0.015;
];

Ksg = governordata(:,1);
K1 = governordata(:,2);
K3 = governordata(:,3);
K5 = governordata(:,4);
K7 = governordata(:,5);
T1 = governordata(:,6);
T2 = governordata(:,7);
T3 = governordata(:,8);
T4 = governordata(:,9);
T5 = governordata(:,10);
T6 = governordata(:,11);
T7 = governordata(:,12);
KAVR = governordata(:,13);
TAVR = governordata(:,14);


GSGInit.Ksg = Ksg;
GSGInit.K1 = K1;
GSGInit.K3 = K3;
GSGInit.K5 = K5;
GSGInit.K7 = K7;
GSGInit.T1 = T1;
GSGInit.T2 = T2;
GSGInit.T3 = T3;
GSGInit.T4 = T4;
GSGInit.T5 = T5;
GSGInit.T6 = T6;
GSGInit.T7 = T7;
GSGInit.KAVR = KAVR;
GSGInit.TAVR = TAVR;


%% 
Pg   = resultsac.gen{:,2}./1000;
Qg   = resultsac.gen{:,3}./1000;
Vg   = resultsac.bus{Syn_Machidx_g,8};
Theta= resultsac.bus{Syn_Machidx_g,9}*pi/180;


Ug=Vg.*exp(1j*Theta);
Sg=Pg+1j*Qg;
IOxy0=conj(Sg./Ug);
EQ=Ug+(1j*Xq).*IOxy0;

Vd = Vg.*sin(angle(EQ)-angle(Ug));
Vq = Vg.*cos(angle(EQ)-angle(Ug));
Iq = Vd./Xq;
Id = abs(IOxy0).*sin(angle(EQ) - angle(IOxy0));
theta0 = angle(EQ) - pi/2; 

Pref = Pg;
% Qref = Qg;

Eq10 = Vq + Id.*Xd1;

Efd0 = Eq10 + Id.*(Xd-Xd1);
Vref0  = Efd0./KAVR+Vg;

GSGInit.Pref = Pref;
GSGInit.Vref0 = Vref0;
GSGInit.theta0 = theta0;
GSGInit.Eq10 = Eq10;
GSGInit.Efd0 = Efd0;
GSGInit.EQ = EQ;

end
