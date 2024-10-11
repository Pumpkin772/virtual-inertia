%% Common parameters
clear
clc

w         =          1;
fB        =          50;
wB        =          2*pi*fB;
W0 = wB;


GFM_in_gidx  =  [3];%The node number of the synchronous motor connected to the large power grid. 
                                    % Please arrange it in the order in the gen matrix.   
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
GFMInit=GFM_init(opdata,GFM_in_gidx);
GFLInit=GFL_init(opdata,GFL_in_gidx);

[Ybus, Yf, Yt] = makeYbus(mpcac);

G_in_gidx = [1,GFL_in_gidx,GFM_in_gidx];



V = opdata{cellstr(string(G_in_gidx)),4};

RG = 0.001;
VO=V(1)+(opdata{1,2}/1000-1i*opdata{1,3}/1000)*RG;

del = opdata{cellstr(string(G_in_gidx)),5}*pi/180;
Vm = V.*cos(del) + 1i*V.*sin(del);
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
[A, B, C, D, YbusPU] = ssNetwMatPower(mpcac, wB);
% 求取状态变量初值
Nu = length(G_in_gidx);
Ne = length(A)/3;
Afull = full(A);
Bfull = full(B);
Adq = Afull(1:2*Ne,1:2*Ne);
Bdq = Bfull(1:2*Ne,1:2*Nu);
udq = [real(Vm);imag(Vm)];
X0=-Adq\Bdq*udq;
X0=[X0;zeros(Ne,1)];
% dx/dt=Ax+Bu u = [Vd(t); Vq(t); V0(t)] (vector 3Nx1).
% y    =Cx+Du y = [Id(t); Iq(t); I0(t)] (vector 3Nx1).
% 更换输入输出
% D^(-1)y=D^(-1)Cx+u 
% u=-D^(-1)Cx + D^(-1)y
% dx/dt=(A-BD^(-1)C)x+BD^(-1)y
A_ = A - B/D*C;
B_ = B/D;
C_ = -D\C;
D_ = inv(D);

units_to_net = [1:3:3*Nu , 2:3:3*Nu , 3:3:3*Nu];
[~,net_to_units] = sort(units_to_net);


%%
AL = [0 -GFLInit.KIL*V(GFL_in_gidx);1 -GFLInit.KPL*V(GFL_in_gidx)];
BL = [-GFLInit.KIL*imag(Vm(GFL_in_gidx))/V(GFL_in_gidx) GFLInit.KIL*real(Vm(GFL_in_gidx))/V(GFL_in_gidx);
    -GFLInit.KPL*imag(Vm(GFL_in_gidx))/V(GFL_in_gidx) GFLInit.KPL*real(Vm(GFL_in_gidx))/V(GFL_in_gidx)];
CL = [0 -GFLInit.Iy0;0 GFLInit.Ix0];

L1 = ACgirdinit.linedata{1, 1}{1,4};
L2 = ACgirdinit.linedata{1, 1}{2,4};
sigma = 1/(1/L1+1/(L2+0.1));
Asys = (eye(2)-sigma/wB*BL*CL)\(AL+sigma*BL*[0 -1;1 0]*CL);
EIG_LINEAR = eig(Asys)

%% 时域模型仿真结果
mdl='OneGFLOneGFM';
open_system(mdl);
io(1)=linio('OneGFLOneGFM/GFL1',1);
io(2)=linio('OneGFLOneGFM/GFM1',1,'output');
linsys1=linearize(mdl,io);
EIG_Complex = eig(linsys1);
EIG_PLL = EIG_Complex(abs(imag(EIG_Complex))>=90&abs(imag(EIG_Complex))<=100)
GFMInit.theta0/2/pi*360