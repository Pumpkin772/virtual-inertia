function [GSGInit] = SG_init(opdata,Syn_Machidx_g)


%%%  Bus No. Xd  Xd1  Xq  Td1 J D
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SGdata=...
[
    1   1.81  0.3   1.78 8 8 20;
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

Pg   =opdata{cellstr(string(Syn_Machidx_g)),2}./1000;
Qg   =opdata{cellstr(string(Syn_Machidx_g)),3}./1000;
Vg   =opdata{cellstr(string(Syn_Machidx_g)),4};
Theta=opdata{cellstr(string(Syn_Machidx_g)),5}*pi/180;

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
