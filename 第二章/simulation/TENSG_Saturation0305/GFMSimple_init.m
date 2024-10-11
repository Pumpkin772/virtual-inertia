function [GFMInit] = GFMSimple_init(resultsac,GFM_in_gidx)

NumGFM = length(GFM_in_gidx);
%%%  Bus No. LG RG CF LF RF J D
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gfmdata=...
[
    30   1e-1   0.0   0.05  0.05 0.01 1 0;
    31   1e-1   0.0   0.05  0.05 0.01 1 0;
    32   1e-1   0.0   0.05  0.05 0.01 1 0;
    33   1e-1   0.0   0.05  0.05 0.01 1 0;
    34   1e-1   0.0   0.05  0.05 0.01 1 0;
    35   1e-1   0.0   0.05  0.05 0.01 1 0;
    36   1e-1   0.0   0.05  0.05 0.01 1 0;
    37   1e-1   0.0   0.05  0.05 0.01 1 0;
    38   1e-1   0.0   0.05  0.05 0.01 1 0;
    39   1e-1   0.0   0.05  0.05 0.01 1 0;
];

LG = gfmdata(1:NumGFM,2);
J  = gfmdata(1:NumGFM,7);
D  = gfmdata(1:NumGFM,8);

GFMInit.LG = LG;
GFMInit.J  = J;
GFMInit.D  = D;

%% 
Pg   = 0;
Qg   = 0;
Vg   = resultsac.bus{GFM_in_gidx,8};
Theta= resultsac.bus{GFM_in_gidx,9}*pi/180;

Ug=Vg.*exp(1j*Theta);
Sg=Pg+1j*Qg;
IOxy0=conj(Sg./Ug);
EQ = Ug + (1j*LG).*IOxy0;
theta0 = angle(EQ) - pi/2; 


GFMInit.Pref = Pg;
GFMInit.Qref = Qg;
GFMInit.EQ   = abs(EQ);
GFMInit.theta0 = theta0;

end
