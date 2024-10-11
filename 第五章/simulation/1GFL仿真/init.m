%% Common parameters
clear
clc

w         =          1;
fB        =          50;                   
wB        =          2*pi*fB;
W0 = wB;


GFM_in_gidx = [1];
GFL_in_gidx = [2];






mpcdc =loadcasedc('P_DCdata');
mpcac =loadcase('P_ACdata');




[resultsac,resultsdc,converged]=powerflowcaculate(mpcac,mpcdc,GFM_in_gidx,GFL_in_gidx);

  

if(converged~=1)
    error('潮流计算不收敛')
end
clear converged

SgridB    =     resultsac.baseMVA;
SdcgridB  =     resultsdc.baseMVAdc;




%% AC network initialization 

fbranchac = resultsac.branch(1,:);  % Set up short line 

[ACgirdinit,opdata]=ACgrid_init(resultsac,resultsdc,fbranchac);

%% Synchronous motor initialization 
GFLInit=GFL_init(opdata,GFL_in_gidx);

[Ybus, Yf, Yt] = makeYbus(mpcac);

V = opdata{cellstr(string(GFL_in_gidx)),4};
del = opdata{cellstr(string(GFL_in_gidx)),5}*pi/180;
Vm = V.*cos(del) + 1i*V.*sin(del);

Vc = 1;

Y = full(Ybus);
Iij = -(Vm - Vc)*Y(2,1);



AL = [0 -GFLInit.KIL*V;1 -GFLInit.KPL*V];
BL = [-GFLInit.KIL*imag(Vm)/V GFLInit.KIL*real(Vm)/V;
    -GFLInit.KPL*imag(Vm)/V GFLInit.KPL*real(Vm)/V];
CL = [0 -GFLInit.Iy0;0 GFLInit.Ix0];

sigma = ACgirdinit.linedata{1, 1}{1,4};
Asys = (eye(2)-sigma/wB*BL*CL)\(AL+sigma*BL*[0 -1;1 0]*CL);
EIG_LINEAR = eig(Asys)

%% 时域模型仿真结果
mdl='OneGFL';
open_system(mdl);
io(1)=linio('OneGFL/GFL1',1);
io(2)=linio('OneGFL/C network',1,'output');
linsys1=linearize(mdl,io);
EIG_Complex = eig(linsys1);
EIG_PLL = EIG_Complex(abs(imag(EIG_Complex))>=90&abs(imag(EIG_Complex))<=500)


