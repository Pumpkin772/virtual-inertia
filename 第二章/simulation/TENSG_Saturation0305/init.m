%% Common parameters
clear
clc

w         =          1;
fB        =          50;                   
wB        =          2*pi*fB;
W0 = wB;


SG_in_gidx     =    30:39;%The node number of the synchronous motor connected to the large power grid. 
                                    % Please arrange it in the order in the gen matrix.   
CSno = SG_in_gidx;% 发电设备

NG = length(SG_in_gidx);

mpcdc =loadcasedc('P_DCdata');
mpcac =loadcase('P_ACdata');

[resultsac,resultsdc,converged]=powerflowcaculate(mpcac,mpcdc,SG_in_gidx);

if(converged~=1)
    error('潮流计算不收敛')
end
clear converged

SgridB    =     resultsac.baseMVA;
SdcgridB  =     resultsdc.baseMVAdc;

% load('twoGFMpara.mat');
load('fourGFMpara.mat');
node = 1:39;

GFM_in_gidx=node(Jes_g.*Ses_g>0.01);
J = Jes_g.*Ses_g;
J = J(GFM_in_gidx);
Satura = Ses_g(GFM_in_gidx);
%% Synchronous motor initialization 
GSGInit=SG_init(resultsac,SG_in_gidx);

GFMInit=GFMSimple_init(resultsac,GFM_in_gidx);

GFMInit.J=J;

% 负荷节点
load_in_gidx=mpcac.bus(mpcac.bus(:,3)~=0|mpcac.bus(:,4)~=0,1);

[Ybus, Yf, Yt] = makeYbus(mpcac);
Vbus = resultsac.bus{:,8}.*exp(1i*resultsac.bus{:,9}*pi/180);

CZno = [];% 恒阻抗负荷
CPno = load_in_gidx';% 恒功率负荷

Vload = Vbus([CZno,CPno]);
Pload = -resultsac.bus{[CZno,CPno],3}/1000;
Qload = -resultsac.bus{[CZno,CPno],4}/1000;

nbus = length(Ybus);


%% 机电暂态
CSno = SG_in_gidx;% 发电设备
CIno = setdiff(1:nbus,unique([CSno,CZno,CPno,GFM_in_gidx]));
U = resultsac.bus{:,8}.*exp(1i*(resultsac.bus{:,9})/180*pi);
I = -(resultsac.bus{:,3}-1i*resultsac.bus{:,4})./conj(U)/1000;
YL = -I(CZno)./U(CZno);
[I_new,Y_new]=kron(Ybus,CZno,CIno,I,YL);
GFMother=setdiff(unique([CSno,CPno,GFM_in_gidx]),[CSno,CPno]);
Y_new = Y_new([CSno,CPno,GFMother],[CSno,CPno,GFMother]);
%% 