function [GFMInit] = GFM_init(opdata,Syn_Machidx_g,SGFM)


%%%  Bus No. LG RG CF LF RF J D
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gfmdata=...
[
    1   0.06   0.0   0.06  0.05 0.001 0.2 10;
];

LG = gfmdata(:,2)./SGFM;
RG = gfmdata(:,3)./SGFM;
CF = gfmdata(:,4).*SGFM;
LF = gfmdata(:,5)./SGFM;
RF = gfmdata(:,6)./SGFM;
J  = gfmdata(:,7).*SGFM;
D  = gfmdata(:,8).*SGFM;

GFMInit.LG = LG;
GFMInit.RG = RG;
GFMInit.CF = CF;
GFMInit.LF = LF;
GFMInit.RF = RF;
GFMInit.J  = J;
GFMInit.D  = D;

%% 
%  KPV KIV KIF KPC KIC KVF TVF
PIcontorldata=...
[
2   10  1 0.3   10  1 0.02;
];

KPV = PIcontorldata(:,1).*SGFM;
KIV = PIcontorldata(:,2).*SGFM;
KIF = PIcontorldata(:,3);
KPC = PIcontorldata(:,4)./SGFM;
KIC = PIcontorldata(:,5)./SGFM;
KVF = PIcontorldata(:,6);
TVF = PIcontorldata(:,7);

GFMInit.KPV = KPV;
GFMInit.KIV = KIV;
GFMInit.KIF = KIF;
GFMInit.KPC = KPC;
GFMInit.KIC = KIC;
GFMInit.KVF = KVF;
GFMInit.TVF = TVF;

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

temp = Idq0 - KIF.*IOdq0 - Vd0.*(1i*CF);
Idref0 = real(temp);
Iqref0 = imag(temp);

GFMInit.Pref = Pref;
GFMInit.Qref = Qref;
GFMInit.Vd0 = Vd0;
GFMInit.theta0 = theta0;
GFMInit.Ix0 = Ix0;
GFMInit.Iy0 = Iy0;
GFMInit.Vx0 = Vx0;
GFMInit.Vy0 = Vy0;
GFMInit.IOx0 = IOx0;
GFMInit.IOy0 = IOy0;
GFMInit.Idref0 = Idref0;
GFMInit.Iqref0 = Iqref0;
GFMInit.Ed0 = Ed0;
GFMInit.Eq0 = Eq0;

end
