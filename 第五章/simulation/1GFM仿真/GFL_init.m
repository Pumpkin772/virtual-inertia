function [GFLInit] = GFL_init(opdata,Syn_Machidx_g)


%%%  Bus No. LG RG CF LF RF J D
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gfmdata=...
[
    2   0.0   0.0   0.0  0.05 0.001 0 0;
    2   0.0   0.0   0.0  0.05 0.001 0 0;
    2   0.0   0.0   0.0  0.05 0.001 0 0;
    2   0.0   0.0   0.0  0.05 0.001 0 0;
    2   0.0   0.0   0.0  0.05 0.001 0 0;
    2   0.0   0.0   0.0  0.05 0.001 0 0;
    2   0.0   0.0   0.0  0.05 0.001 0 0;
    2   0.0   0.0   0.0  0.05 0.001 0 0;
    2   0.0   0.0   0.0  0.05 0.001 0 0;
];

LG = gfmdata(:,2);
RG = gfmdata(:,3);
CF = gfmdata(:,4);
LF = gfmdata(:,5);
RF = gfmdata(:,6);
JGFL  = gfmdata(:,7);
DGFL  = gfmdata(:,8);

GFLInit.LG = LG;
GFLInit.RG = RG;
GFLInit.CF = CF;
GFLInit.LF = LF;
GFLInit.RF = RF;
GFLInit.J  = JGFL;
GFLInit.D  = DGFL;

%% 
%  KPL KIL TFL KPP KIP TPF KPC KIC KVF TVF
PIcontorldata=...
[
8 9500  0.01 0.1 5  0.05 2   10  1 0.001;
8 9500  0.01 0.1 5  0.05 2   10  1 0.001;
8 9500  0.01 0.1 5  0.05 2   10  1 0.001;
8 9500  0.01 0.1 5  0.05 2   10  1 0.001;
8 9500  0.01 0.1 5  0.05 2   10  1 0.001;
8 9500  0.01 0.1 5  0.05 2   10  1 0.001;
8 9500  0.01 0.1 5  0.05 2   10  1 0.001;
8 9500  0.01 0.1 5  0.05 2   10  1 0.001;
8 9500  0.01 0.1 5  0.05 2   10  1 0.001;
];

KPL = PIcontorldata(:,1);
KIL = PIcontorldata(:,2);
TFL = PIcontorldata(:,3);
KPP = PIcontorldata(:,4);
KIP = PIcontorldata(:,5);
TPF = PIcontorldata(:,6);
KPC = PIcontorldata(:,7);
KIC = PIcontorldata(:,8);
KVF = PIcontorldata(:,9);
TVF = PIcontorldata(:,10);

GFLInit.KPL = KPL;
GFLInit.KIL = KIL;
GFLInit.TFL = TFL;
GFLInit.KPP = KPP;
GFLInit.KIP = KIP;
GFLInit.TPF = TPF;
GFLInit.KPC = KPC;
GFLInit.KIC = KIC;
GFLInit.KVF = KVF;
GFLInit.TVF = TVF;

%% 

Pg   =opdata{cellstr(string(Syn_Machidx_g)),2}./1000;
Qg   =opdata{cellstr(string(Syn_Machidx_g)),3}./1000;
Vg   =opdata{cellstr(string(Syn_Machidx_g)),4};
Theta=opdata{cellstr(string(Syn_Machidx_g)),5}*pi/180;

Ug=Vg.*exp(1j*Theta);
Sg=Pg+1j*Qg;
IOxy0=conj(Sg./Ug);
Vxy0=Ug+(RG+1j*LG).*IOxy0;
Ixy0  = IOxy0 + Vxy0.*(1i*CF);
Exy0  = Vxy0 + Ixy0.*(RF+1i*LF);

IOx0 = real(IOxy0);
IOy0 = imag(IOxy0);
Vx0 = real(Vxy0);
Vy0 = imag(Vxy0);
Ix0 = real(Ixy0);
Iy0 = imag(Ixy0);
Pref = real(Vxy0.*conj(IOxy0));
Qref = imag(Vxy0.*conj(IOxy0));
Vd0  = abs(Vxy0);
theta0 = angle(Vxy0);

Edq0 = abs(Exy0).*exp(1i*(angle(Exy0)-theta0));
Idq0 = abs(Ixy0).*exp(1i*(angle(Ixy0)-theta0));
IOdq0 = abs(IOxy0).*exp(1i*(angle(IOxy0)-theta0));

temp = Edq0 - KVF.*Vd0 - Idq0.*(1i*LF);
Ed0 = real(temp);
Eq0 = imag(temp);

Idref0 = real(Idq0);
Iqref0 = imag(Idq0);

GFLInit.Pref = Pref;
GFLInit.Qref = Qref;
GFLInit.Vd0 = Vd0;
GFLInit.theta0 = theta0;
GFLInit.Ix0 = Ix0;
GFLInit.Iy0 = Iy0;
GFLInit.Vx0 = Vx0;
GFLInit.Vy0 = Vy0;
GFLInit.IOx0 = IOx0;
GFLInit.IOy0 = IOy0;
GFLInit.Idref0 = Idref0;
GFLInit.Iqref0 = Iqref0;
GFLInit.Ed0 = Ed0;
GFLInit.Eq0 = Eq0;

end
