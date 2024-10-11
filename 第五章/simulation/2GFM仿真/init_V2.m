%% Common parameters
clear
clc

w         =          1;
fB        =          50;                   
wB        =          2*pi*fB;
W0 = wB;


% GFM_in_gidx  =  [3];%The node number of the synchronous motor connected to the large power grid. 
                                    % Please arrange it in the order in the gen matrix.   
% GFL_in_gidx = [1 2 3 5 7:11]; 
GFL_in_gidx = [];
SG_in_gidx  = [39];
GFM_in_gidx = [4 6];







mpcdc =loadcasedc('P_DCdata');
mpcac =loadcase('P_ACdata');

% Pref = [100;500;1000;1000;700;1000;600;1000;1000];

% Pref = [100;500;400;300;700;1000;600;1000;200];
% Pref = 100*ones(9,1);
% Pref = [800;900;700;800;1000;700;500;1000;800];
% Pref = 1000*ones(9,1);
% mpcac.gen(GFL_in_gidx,2)=Pref;


[resultsac,resultsdc,converged]=powerflowcaculate(mpcac,mpcdc,SG_in_gidx,GFM_in_gidx,GFL_in_gidx);

  

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
SGFM = [1;1];
GFMInit=GFM_init(opdata,GFM_in_gidx,SGFM);
% GFLInit=GFL_init(opdata,GFL_in_gidx);

% GFLInit.KPL=[8;8;8;8;9;8;9;8;9];
% GFLInit.KIL=[9500;9500;9500;9500;9000;9500;9000;9500;9000];

[Ybus, Yf, Yt] = makeYbus(mpcac);

G_in_gidx = [SG_in_gidx,GFL_in_gidx,GFM_in_gidx];

G_in_gidx = sort(G_in_gidx);

Vm = resultsac.bus{:,8}.*exp(1i*resultsac.bus{:,9}*pi/180);
V  = Vm;

RG = 1e-6;
VO=Vm(SG_in_gidx)+(opdata{num2str(SG_in_gidx),2}/1000-1i*opdata{num2str(SG_in_gidx),3}/1000)*RG;



fb = resultsac.branch{:,1};                % From bus number...
tb = resultsac.branch{:,2};                % To bus number...
nl = length(fb);                   % No. of Branches..
Pl = resultsac.bus{:,3};                 % PLi..
Ql = resultsac.bus{:,4};                 % QLi..

nbus = length(Ybus);
Iij = zeros(nbus,nbus);
Si = zeros(nbus,1);

%Line Current Flows..
for m = 1:nl
    p = fb(m); q = tb(m);
    Iij(p,q) = -(Vm(p) - Vm(q))*Ybus(p,q); % Y(m,n) = -y(m,n)..
    Iij(q,p) = -Iij(p,q);
end

% Compute dq0 state-space model
[A1, B1, C1, D1, YbusPU] = ssNetwMatPower(mpcac, wB);
subset = G_in_gidx;

[A,B,C,D] = elimBus(A1,B1,C1,D1,subset); 
% 求取状态变量初值
Nu = length(G_in_gidx);
Ne = length(A)/3;
Afull = full(A);
Bfull = full(B);
Adq = Afull(1:2*Ne,1:2*Ne);
Bdq = Bfull(1:2*Ne,1:2*Nu);
udq = [real(Vm(G_in_gidx));imag(Vm(G_in_gidx))];
X0=-Adq\Bdq*udq;
% X0=[X0;zeros(Ne,1)];
% dx/dt=Ax+Bu u = [Vd(t); Vq(t); V0(t)] (vector 3Nx1).
% y    =Cx+Du y = [Id(t); Iq(t); I0(t)] (vector 3Nx1).
% 更换输入输出
% D^(-1)y=D^(-1)Cx+u 
% u=-D^(-1)Cx + D^(-1)y
% dx/dt=(A-BD^(-1)C)x+BD^(-1)y
A_ = A(1:2*Ne,1:2*Ne) - B(1:2*Ne,1:2*Nu)/D(1:2*Nu,1:2*Nu)*C(1:2*Nu,1:2*Ne);
B_ = B(1:2*Ne,1:2*Nu)/D(1:2*Nu,1:2*Nu);
C_ = -D(1:2*Nu,1:2*Nu)\C(1:2*Nu,1:2*Ne);
D_ = inv(D(1:2*Nu,1:2*Nu));

units_to_net = [1:2:2*Nu , 2:2:2*Nu];
[~,net_to_units] = sort(units_to_net);


%%

% 
% 


% 
% L1 = ACgirdinit.linedata{1, 1}{1,4};
% L2 = ACgirdinit.linedata{1, 1}{2,4};
% sigma = 1/(1/L1+1/(L2+0.1));

% 
% Q = zeros(39,39);
% Q(1,11) = -92.08;
% Q(2,15) = -66.67;
% Q(3,19) = -83.33;
% Q(4,28) = -117.37;
% Q(5,29) = -92.59;
% Q(6,31) = -116.55;
% Q(7,34) = -71.84;
% Q(8,38) = -106.84;
% Q(9,18) = -66.67;
% Q(10,11) = -40.55;
% Q(11,12) = -110.38;
% Q(12,13) = -78.25;
% Q(13,14) = -130.21;
% Q(14,15) = -641.03;
% Q(15,16) = -181.16;
% Q(16,17) = -362.32;
% Q(17,18) = -45.91;
% Q(15,20) = -203.25;
% Q(19,20) = -387.60;
% Q(20,21) = -38.31;
% Q(19,22) = -387.60;
% Q(21,22) = -38.31;
% Q(13,23) = -129.20;
% Q(22,23) = -165.02;
% Q(23,24) = -76.80;
% Q(24,25) = -177.31;
% Q(25,26) = -187.27;
% Q(12,27) = -125.31;
% Q(26,27) = -203.25;
% Q(25,28) = -85.47;
% Q(28,29) = -120.77;
% Q(25,30) = -123.46;
% Q(30,31) = -119.05;
% Q(31,32) = -173.61;
% Q(25,33) = -282.49;
% Q(32,33) = -47.62;
% Q(11,34) = -193.80;
% Q(34,35) = -52.60;
% Q(26,36) = -96.34;
% Q(35,36) = -113.38;
% Q(35,37) = -35.16;
% Q(35,38) = -26.67;
% Q(37,38) = -110.38;
% Q(32,39) = -61.27;
% Q(9,10) = -66.67;
% Q(14,17) = -148.81;
% Q = Q + Q';
% Q = -diag(sum(Q))+Q;
% Q(39,:) = [];
% Q(:,39) = [];
% 
% 
% 
% L_ = setdiff(1:38,GFL_in_gidx);
% Qred = Q(GFL_in_gidx,GFL_in_gidx) -  Q(GFL_in_gidx,L_)*Q(L_,L_)^(-1)*Q(L_,GFL_in_gidx);
% [VV,DD,WW] = eig(Qred);
% [lamb1,inx] = min(diag(DD));
% pp=VV(:,3).*VV(:,3);
% AL = 0;
% BL = 0;
% CL = 0;
% for i = 1:9
%     AL = AL + pp(i)*[0 -GFLInit.KIL(i)*abs(V(GFL_in_gidx(i)));1 -GFLInit.KPL(i)*abs(V(GFL_in_gidx(i)))];
%     BL = BL + pp(i)*[-GFLInit.KIL(i)*imag(V(GFL_in_gidx(i)))/abs(V(GFL_in_gidx(i))) GFLInit.KIL(i)*real(V(GFL_in_gidx(i)))/abs(V(GFL_in_gidx(i)));
%         -GFLInit.KPL(i)*imag(V(GFL_in_gidx(i)))/abs(V(GFL_in_gidx(i))) GFLInit.KPL(i)*real(V(GFL_in_gidx(i)))/abs(V(GFL_in_gidx(i)))];
%     CL = CL + pp(i)*[0 -GFLInit.Iy0(i);0 GFLInit.Ix0(i)];
% end
% sigma = 1/lamb1;
% Asys = (eye(2)-sigma/wB*BL*CL)\(AL+sigma*BL*[0 -1;1 0]*CL);
% EIG_LINEAR = eig(Asys)
% 
% %% 
% [VV2,DD2,WW2] = eig(diag((abs(V(GFL_in_gidx)).^2)./(Pref/1000))*Qred);
% [lamb1_2,inx_2] = min(diag(DD2));
% pp2=VV(:,inx_2).*VV(:,inx_2);
% Ueq = 1 + pp2'*(abs(V(GFL_in_gidx))-1);
% syms s
% HPLL = (GFLInit.KPL(1)+GFLInit.KIL(1)/s)/s;
% Yxygfl = [0 0;0 -Ueq*HPLL/(1+Ueq*HPLL)];
% tao = 0.0;
% gama = 1/((s+tao)^2/W0+W0)*[s+tao W0;-W0 s+tao];
% double(solve(det(Yxygfl+kron(lamb1_2,gama)),s))
% 
% HPLL = (GFLInit.KPL(2)+GFLInit.KIL(2)/s)/s;
% Yxygfl = [0 0;0 -Ueq*HPLL/(1+Ueq*HPLL)];
% tao = 0.0;
% gama = 1/((s+tao)^2/W0+W0)*[s+tao W0;-W0 s+tao];
% double(solve(det(Yxygfl+kron(lamb1_2,gama)),s))
    
%% 时域模型仿真结果
mdl='NineGFL_V2';
open_system(mdl);
io(1)=linio('NineGFL_V2/GFM1',1);
io(2)=linio('NineGFL_V2/transmission network',1,'output');
linsys1=linearize(mdl,io);
[EIG_Complex] = eig(linsys1);
EIG_PLL = EIG_Complex(abs(imag(EIG_Complex))>=90&abs(imag(EIG_Complex))<=100)