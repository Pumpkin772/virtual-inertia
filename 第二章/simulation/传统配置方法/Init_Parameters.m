%% begin the program
clear;clc;
%% Parameters Initialize
%% Case Set
% load case
% mpc=loadcase('P_ACdata');

%% Base Capacity
% baseMVA=mpc.baseMVA;
% bus=mpc.bus;
% branch=mpc.branch;
% DC power flow sensitivity matrix
% S_b=makePTDF(baseMVA,bus,branch);

%% 不同位置配置GFM储能后的网架方程
% 假设节点7、节点9为待配置网络节点，还包括节点1,2,3,4发电侧
mpcdc =loadcasedc('P_DCdata');
mpcac =loadcase('P_ACdata');

% mpcdc=loadcasedc('P_DCdata');
% 基态潮流
[resultsac,resultsdc,converged] = powerflowcaculate(mpcac,mpcdc);
U = resultsac.bus{:,8}.*exp(1i*(resultsac.bus{:,9})/180*pi);
I = -(resultsac.bus{:,3}-1i*resultsac.bus{:,4})./conj(U)/1000;
[Ybus, Yf, Yt] = makeYbus(mpcac);
Nbus = length(Ybus);
baseMVA = 1000; % 基准容量，如果要做标幺值有名值换算就有用

NodeSG = [30:39];
GSGInit=SG_init(resultsac,NodeSG);

%% 节点6和节点10是否具备构网储能，会影响转移矩阵
% 分别计算不同组合情况下的转移矩阵

Pload = -resultsac.bus{:,3}/1000; % 基准值为1000MW，设置四个扰动点
deltaPD = 0.2*Pload(Pload~=0);
fau = find(Pload~=0);% 节点1、2、3、4功率扰动
% fau = 3;
NodeESCandidate = [1:39]+10+39;
% NodeESCandidate = [21:22] + 10+39;
NodeESCandidateCount = length(NodeESCandidate);

NodeSG_new = NodeSG+10;
NodeS = NodeSG_new;
NodeS = unique([NodeESCandidate,NodeS]);
GenCount = length(NodeS);

NESCan   = length(NodeESCandidate);
Nbus_new = Nbus + 10 + NESCan;
Ybus_ESCan = zeros(Nbus_new,Nbus_new);
Ybus_ESCan(1:Nbus,1:Nbus) = Ybus;


Ybus_ESCan(NodeSG,NodeSG)=Ybus_ESCan(NodeSG,NodeSG)+diag(1./(1i*GSGInit.Xq));
Ybus_ESCan(NodeSG,NodeSG_new)=-diag(1./(1i*GSGInit.Xq));
Ybus_ESCan(NodeSG_new,NodeSG)=-diag(1./(1i*GSGInit.Xq));
Ybus_ESCan(NodeSG_new,NodeSG_new)=diag(1./(1i*GSGInit.Xq));

Ybus_ESCan(1:NESCan,Nbus+11:Nbus_new) = -1/(0.1*1i)*eye(length(NodeESCandidate));
Ybus_ESCan(Nbus+11:Nbus_new,1:NESCan) = -1/(0.1*1i)*eye(length(NodeESCandidate));
Ybus_ESCan(Nbus+11:Nbus_new,Nbus+11:Nbus_new) = 1/(0.1*1i)*eye(length(NodeESCandidate));
Ybus_ESCan(1:NESCan,1:NESCan) = Ybus_ESCan(1:NESCan,1:NESCan) + 1/(0.1*1i)*eye(length(NodeESCandidate));


I = [I;zeros(10,1);zeros(NESCan,1)];
U = [U;GSGInit.EQ;U(1:NESCan)];
del_S = angle(U);
L0=zeros(length(U),length(U));
G = real(Ybus_ESCan);          	% Conductance matrix..
B = imag(Ybus_ESCan);         	% Susceptance matrix..
for m = 1:length(U)
    for n = 1:length(U)
        if m ~= n
            L0(m,n) = -abs(U(m))*abs(U(n))*(-G(m,n)*sin(del_S(m)-del_S(n))+B(m,n)*cos(del_S(m)-del_S(n)));
        end
    end
    L0(m,m) = -sum(L0(m,:));
end

%% 节点削减
% MaxCombination = 1+ NodeESCandidateCount+nchoosek(NodeESCandidateCount,2); % 选择两个位置配置储能
% Combination = cell(MaxCombination,1);
% Combination{1} = [];
% for i = 1:NodeESCandidateCount
%     Combination{i+1} = NodeESCandidate(i);
% end
% Comb = nchoosek(NodeESCandidate,2);
% for i = 1:nchoosek(NodeESCandidateCount,2)
%     Combination{i+1+NodeESCandidateCount} = Comb(i,:);
% end
NodeESCandidate = [19:22] + 10+39;
NodeS = [NodeSG_new,NodeESCandidate];
NodeESCandidateCount = length(NodeESCandidate);
MaxCombination = 2^NodeESCandidateCount;
ind=ff2n(NodeESCandidateCount)>0;
Combination = arrayfun(@(i)NodeESCandidate(:,ind(i,:)),1:MaxCombination,'UniformOutput',false);
Y = cell(1,MaxCombination);
L = cell(1,MaxCombination);

for i = 1:MaxCombination
NodeCS = [NodeSG_new,Combination{i}];
NodeCZ = [];% 恒阻抗负荷
NodeCP = fau';
NodeCI = setdiff(1:Nbus_new,unique([NodeCS,NodeCZ,NodeCP]));
YL = -I(NodeCZ)./U(NodeCZ);
[I_new,Y_new]=kron_user(Ybus_ESCan,NodeCZ,NodeCI,I,YL);
Node_re=unique([NodeCP,NodeS]);
NodeCount=length(Node_re);
Y{i}=Y_new(Node_re,Node_re);
L{i}=zeros(NodeCount,NodeCount);
U_S = abs(U(Node_re));
del_S = angle(U(Node_re));
G = real(Y{i});          	% Conductance matrix..
B = imag(Y{i});         	% Susceptance matrix..
for m = 1:NodeCount
    for n = 1:NodeCount
        if m ~= n
            L{i}(m,n) = -U_S(m)*U_S(n)*(-G(m,n)*sin(del_S(m)-del_S(n))+B(m,n)*cos(del_S(m)-del_S(n)));
        end
    end
    L{i}(m,m) = -sum(L{i}(m,:));
end
end
% save('ParaInit')
