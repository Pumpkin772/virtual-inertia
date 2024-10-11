clear;
clc;
tic
Init_Parameters;
J = [4.5 11.6 11.6 4.5 11.6 1.9 1.7 11.6 11.6+1 11.6];

% 储能惯性系数约束
Jemin = 4;
Jemax = 12;


% 频率变化率限制
RoCoFlim = 0.5;
FreqBase = 50;
DeltFlim = 0.5;

% 储能参数
Cp = 300000;% 功率容量成本系数（$/MW）
Ce = 150000;% 能量容量成本系数（$/MWh）
MaxESCount = 2;% 最大储能数量

num_fau = length(deltaPD);          % 扰动数量
% num_fau = 1;
NodeSG  = setdiff(NodeS,NodeESCandidate);
SGCount = length(NodeSG);




%% 定义决策变量   sdpvar定义连续变量，binvar定义整数变量；待求解
Ses_g=sdpvar(NodeESCandidateCount,1); % 网侧储能功率容量，单位pu，1pu=1000MW
Ees_g=sdpvar(NodeESCandidateCount,1); % 网侧储能能量容量，单位MWh
Jes_g=sdpvar(NodeESCandidateCount,1); % 网侧储能惯性系数

BES = binvar(NodeESCandidateCount,1);% 储能数量
SL = binvar(MaxCombination,1); % 选择网络的0-1变量
Ves=sdpvar(1,1); % 选择网络矩阵形式的中间变量，代表待选节点组合形式


LPD=sdpvar(length(NodeCS),length(deltaPD)); % L12L22^(-1)
LPdelta=sdpvar(length(NodeCS),length(NodeCS),'full'); % L11-L12L22^(-1)L21


% 惯性系数与网侧储能功率容量乘积二进制扩展
K = 20;
ves_gk = binvar(NodeESCandidateCount,K,'full');
gamak =  sdpvar(NodeESCandidateCount,K,'full');
delta_Je = (Jemax-Jemin)/(2^K);




Sa_c=cell(num_fau,1);% 不同扰动下的网侧待选节点储能是否饱和0-1变量
Del_c=cell(num_fau,1);% 不同扰动下的网侧待选节点相角变化
% Del=cell(num_fau,1);% 不同扰动下的所有节点的相角变化
Pes_g_c=cell(num_fau,1);% 不同扰动下的网侧储能输出功率
Psg=cell(num_fau,1);% 不同扰动下的发电机输出功率
RoCoF_sg=cell(num_fau,1);% 不同扰动下的各发电机节点频率变化率


Ze_c=cell(num_fau,1);% 不同扰动下的惯性中心点频率变化率约束处理中间0-1变量
Ye_c=cell(num_fau,1);% 不同扰动下的惯性中心点频率变化率约束处理储能承担的功率中间0-1变量
Ee_c=cell(num_fau,1);% 能量约束

for Ind_fau=1:num_fau
    Sa_c{Ind_fau,1}=binvar(NodeESCandidateCount,1); % 网侧待选节点储能是否饱和0-1变量
    Del_c{Ind_fau,1}=sdpvar(NodeESCandidateCount,1); % 网侧待选节点相角变化 
%     Del{Ind_fau,1}=sdpvar(NodeCount,1);% 所有节点的相角变化
    Pes_g_c{Ind_fau,1}=sdpvar(NodeESCandidateCount,1); % 网侧储能输出功率
    Psg{Ind_fau,1}=sdpvar(SGCount,1);
    RoCoF_sg{Ind_fau,1}=sdpvar(SGCount,1);
    
    Ze_c{Ind_fau,1}=sdpvar(NodeESCandidateCount,1);
    Ye_c{Ind_fau,1}=sdpvar(NodeESCandidateCount,1);
    Ee_c{Ind_fau,1}=sdpvar(NodeESCandidateCount,1);
end


%% 约束条件——系统
Constraints = [];
Constraints = [Constraints;sum(BES)<= MaxESCount];% 储能数量约束
Constraints = [Constraints;0<=Ses_g<=BES];
% 待选节点组合If-then logic
Uweight = [];
for i=1:NodeESCandidateCount
    Uweight=[2^(i-1),Uweight];
end
Constraints = [Constraints;Uweight*BES+1==Ves];      % 组合形式的二进制数
Constraints = [Constraints;sum(SL)==1]; % 选择网络只能有一个
Constraints = [Constraints;(1:MaxCombination)*SL==Ves];
[~,Distlocation]=ismember(fau',Node_re);
[~,Inerlocation]=ismember(setdiff(Node_re,fau'),Node_re);
for i = 1:MaxCombination
    Constraints = [Constraints;-10^6*(1-SL(i,1))<=LPD-L{i}(Inerlocation,Distlocation)/L{i}(Distlocation,Distlocation)<=10^6*(1-SL(i,1))];
    Constraints = [Constraints;-10^6*(1-SL(i,1))<=LPdelta-(L{i}(Inerlocation,Inerlocation)-L{i}(Inerlocation,Distlocation)/L{i}(Distlocation,Distlocation)*L{i}(Distlocation,Inerlocation))<=10^6*(1-SL(i,1))];
end
% Constraints = [Constraints;LPD==L{4}(Inerlocation,Distlocation)/L{4}(Distlocation,Distlocation)];
% Constraints = [Constraints;LPdelta==L{4}(Inerlocation,Inerlocation)-L{4}(Inerlocation,Distlocation)/L{4}(Distlocation,Distlocation)*L{4}(Distlocation,Inerlocation)];

% Constraints = [Constraints;0<=Ses_s];
Constraints=[Constraints,Jes_g==Jemin+delta_Je*(ves_gk*(2.^[[1:K]-1])')];% 惯性系数，惯性系数乘以容量表达式为Ses_g*Jemin+delta_Je*(gamak*(2.^[[1:K]-1])')
Constraints=[Constraints,-100*ves_gk<=gamak<=100*ves_gk];% gamak=ves_gk*Ses_g，0-1变量乘以连续变量
for k=1:K
    Constraints=[Constraints,Ses_g-100*(1-ves_gk(:,k))<=gamak(:,k)<=Ses_g+100*(1-ves_gk(:,k))];
end
for Ind_fau = 1:num_fau
    [~,Inerlocation]=ismember(setdiff(NodeS,fau(Ind_fau)),Node_re);
    [~,Distlocation]=ismember(fau(Ind_fau),Node_re);
    [~,Distotherlocation]=ismember(setdiff(NodeCP,fau(Ind_fau)),Node_re);
    [~,ESCanlocation]=ismember(setdiff(NodeESCandidate,fau(Ind_fau)),Node_re);
    [~,ESInterCanlocation]=ismember(setdiff(NodeESCandidate,fau(Ind_fau)),NodeESCandidate);
    [~,Sglocation]=ismember(setdiff(NodeSG,fau(Ind_fau)),Node_re);
    [~,SglocatSou]=ismember(setdiff(NodeSG,fau(Ind_fau)),setdiff(NodeS,NodeESCandidate));
    [~,ESlocatSou]=ismember(setdiff(NodeESCandidate,fau(Ind_fau)),NodeS);
    [~,DistlocatSou]=ismember(fau(Ind_fau),setdiff(NodeS,NodeESCandidate));
    [~,DistlocatES]=ismember(fau(Ind_fau),NodeESCandidate);

    % 扰动节点相角变化
%     Constraints = [Constraints;deltaPD(Ind_fau)==L(Distlocation,:)*Del{Ind_fau,1}];% 扰动节点功率对相角约束
%     Constraints = [Constraints;0==L(Distotherlocation,:)*Del{Ind_fau,1}];% 其他恒功率负载节点
%     Constraints = [Constraints;Del{Ind_fau,1}(Sglocation,1)==0];% 发电机节点相角不变
%     Constraints = [Constraints;Del{Ind_fau,1}(ESCanlocation,1)== Del_c{Ind_fau,1}];% 储能待配置节点
%     Constraints = [Constraints;Del{Ind_fau,1}(Distlocation,1)<=Del{Ind_fau,1}];% 为了保证其他节点功率输出并不任意，扰动节点的相位变化是最大的L22-1的性质
%     Constraints = [Constraints;-2*pi<=Del{Ind_fau,1}<=0*pi];% 这个也是L22-1的性质
    
    deltaP = zeros(num_fau,1);
    deltaP(Ind_fau)=deltaPD(Ind_fau);
    % 储能输出功率约束
    Constraints = [Constraints;0<=Pes_g_c{Ind_fau,1}<=Ses_g]; % 储能输出功率限制
    Constraints = [Constraints;-10*(1-Sa_c{Ind_fau,1})<=Ses_g-Pes_g_c{Ind_fau,1}<=10*(1-Sa_c{Ind_fau,1})]; % 大M法，储能没有达到饱和时，相角为0；达到饱和时，功率为最大输出功率
    Constraints = [Constraints;-2*pi*(Sa_c{Ind_fau,1})*0<=Del_c{Ind_fau,1}...
        <=0*pi*(Sa_c{Ind_fau,1})]; % 大M法，储能没有达到饱和时，相角为0；达到饱和时，功率为最大输出功率
    
    Constraints = [Constraints;Pes_g_c{Ind_fau,1}==LPdelta(ESlocatSou,ESlocatSou)*Del_c{Ind_fau,1}+...
        LPD(ESlocatSou,:)*deltaP];% 储能输出功率和网络之间的关系
    % 发电机节点频率变化率约束
%     Constraints = [Constraints;Psg{Ind_fau,1}==L(Sglocation,:)*Del{Ind_fau,1}];
    Constraints = [Constraints;Psg{Ind_fau,1}==LPdelta(SglocatSou,ESlocatSou)*Del_c{Ind_fau,1}+LPD(SglocatSou,:)*deltaP];
    Constraints = [Constraints;RoCoF_sg{Ind_fau,1}==(-Psg{Ind_fau,1})./J'*FreqBase];% 发电机输出功率除以惯性
    Constraints = [Constraints;-RoCoFlim<=RoCoF_sg{Ind_fau,1}<=RoCoFlim];

        
    % 附加的惯性中心点频率变化率约束
    Constraints = [Constraints;-1000*(1-Sa_c{Ind_fau,1})<=Ze_c{Ind_fau,1}<=1000*(1-Sa_c{Ind_fau,1})];% 储能饱和则无法体现出其惯性，但是相应的扰动会减少
    Constraints = [Constraints;(Pes_g_c{Ind_fau,1}*FreqBase-RoCoFlim*(Ses_g*Jemin+delta_Je*(gamak*(2.^[[1:K]-1])')))-1000*(Sa_c{Ind_fau,1})<=Ze_c{Ind_fau,1}...
        <=(Pes_g_c{Ind_fau,1}*FreqBase-RoCoFlim*(Ses_g*Jemin+delta_Je*(gamak*(2.^[[1:K]-1])')))+1000*(Sa_c{Ind_fau,1})];% 储能饱和则无法体现出其惯性，但是相应的扰动会减少
    Constraints = [Constraints;sum(Psg{Ind_fau,1}*FreqBase-RoCoFlim*J')+sum(Ze_c{Ind_fau,1})<=0];% 惯性中心点频率变化率
    
    Constraints = [Constraints;RoCoFlim/FreqBase*(Ses_g*Jemin+delta_Je*(gamak*(2.^[[1:K]-1])'))-1000*(Sa_c{Ind_fau,1})<=Ye_c{Ind_fau,1}...
        <=RoCoFlim/FreqBase*(Ses_g*Jemin+delta_Je*(gamak*(2.^[[1:K]-1])'))+1000*(Sa_c{Ind_fau,1})]; % 每台储能惯性支撑的功率大小
    Constraints = [Constraints;-1000*(1-Sa_c{Ind_fau,1})<=Ye_c{Ind_fau,1}<=1000*(1-Sa_c{Ind_fau,1})];
    Constraints = [Constraints;Ye_c{Ind_fau,1}<=Ses_g];
    
    % 储能节点频率变化率约束《正常来说没有这个要求，但是如果惯性较小，靠近该处的机组在后续过程中依然会》
    Constraints = [Constraints;Ze_c{Ind_fau,1}+Ye_c{Ind_fau,1}*FreqBase<=RoCoFlim*(Ses_g*Jemin+delta_Je*(gamak*(2.^[[1:K]-1])'))];% Ze_c+Ye_c*f0=(1-s)*Pesg*f0
        % 储能调频能量需求
    Constraints = [Constraints;Ees_g-DeltFlim/RoCoFlim*Ye_c{Ind_fau,1}*1000>=0];
    Constraints = [Constraints;-100*Sa_c{Ind_fau,1}<=Ee_c{Ind_fau,1}<=100*Sa_c{Ind_fau,1}];
    Constraints = [Constraints;Ses_g-100*(1-Sa_c{Ind_fau,1})<=Ee_c{Ind_fau,1}<=100*(1-Sa_c{Ind_fau,1})+Ses_g];
    Constraints = [Constraints;Ees_g>=Ee_c{Ind_fau,1}*5/60*1000];
end
% 频率变化率约束
% Constraints = [Constraints;Ses_s(SglocatSou,1)==0];
%% 目标函数及其求解

Objective = sum(Ses_g)*1000*Cp+sum(Ees_g)*Ce;
options=sdpsettings('solve','gurobi');
options=sdpsettings(options,'verbose',2);
result=solvesdp(Constraints,Objective,options);
toc
%% 变量转换
Ses_g=double(Ses_g);
Jes_g=double(Jes_g);
Ees_g=double(Ees_g);
BES=double(BES);
Objective=double(Objective);

RoCoF_cen = zeros(num_fau,1);
RoCoF_total = cell(num_fau,1);
Psg_total = cell(num_fau,1);
Pes_g_total = cell(num_fau,1);

for Ind_fau=1:num_fau
    Sa_c{Ind_fau,1}=value(Sa_c{Ind_fau,1}); 
    Del_c{Ind_fau,1}=value(Del_c{Ind_fau,1}); 
    Del{Ind_fau,1}=value(Del{Ind_fau,1}); 
    Pes_g_c{Ind_fau,1}=value(Pes_g_c{Ind_fau,1});
    RoCoF_sg{Ind_fau,1}=value(RoCoF_sg{Ind_fau,1});
    Psg{Ind_fau,1}=value(Psg{Ind_fau,1});
    
    Ze_c{Ind_fau,1}=value(Ze_c{Ind_fau,1});
    Ye_c{Ind_fau,1}=value(Ye_c{Ind_fau,1});
    
    RoCoF_total{Ind_fau,1}=zeros(1,Nbus_new);
    Psg_total{Ind_fau,1}=Psg{Ind_fau,1}';
    
    
    [~,SglocatSou]=ismember(setdiff(NodeSG,fau(Ind_fau)),setdiff(NodeS,NodeESCandidate));
    RoCoF_total{Ind_fau,1}(NodeSG)=RoCoF_sg{Ind_fau,1};
    Pes_g_total{Ind_fau,1}=Pes_g_c{Ind_fau,1}';
    
    % 求解各节点频率变化率
    [~,InertiaES_location] = ismember(NodeESCandidate(Sa_c{Ind_fau,1}<1),1:Nbus_new);
    uninertia_node = [NodeCP,NodeESCandidate(Sa_c{Ind_fau,1}==1)];
    inertia_node = setdiff(Node_re,uninertia_node);
    uninertia_node = setdiff(1:Nbus_new,inertia_node);
%     [~,inertia_node_location]=ismember(inertia_node,Node_re);
%     [~,uninertia_node_location]=ismember(uninertia_node,Node_re);
    
    % 储能节点频率变化率
    RoCoF_BESS=-Pes_g_c{Ind_fau,1}.*(1-Sa_c{Ind_fau,1})./(Jes_g.*Ses_g)*FreqBase;
    RoCoF_BESS=RoCoF_BESS(Sa_c{Ind_fau,1}<1);
    RoCoF_total{Ind_fau,1}(InertiaES_location)=RoCoF_BESS;
    % 非惯性节点频率变化率
    RoCoF_uninertia=-inv(L0(uninertia_node,uninertia_node))*L0(uninertia_node,inertia_node)*[RoCoF_sg{Ind_fau,1};RoCoF_BESS];
    RoCoF_total{Ind_fau,1}(uninertia_node)=RoCoF_uninertia;
    
    % 惯性中心点频率变化率
    RoCoF_cen(Ind_fau) = (sum(Psg{Ind_fau,1})+sum(Pes_g_c{Ind_fau,1}.*(1-Sa_c{Ind_fau,1})))/(sum(J)+sum(Jes_g.*Ses_g.*(1-Sa_c{Ind_fau,1})))*FreqBase;
end
% save('TENSG_Saturation0305/oneGFMconfig');
%% 
% 绘制频率变化率
figure(1)
bar3(-cell2mat(RoCoF_total),0.6);
colormap(slanCM(100))
hXLabel = xlabel('机组');
hYLabel = ylabel('场景');
hZLabel = zlabel('频率变化率');

hTitle=title('不同扰动场景下各机组频率变化率');

% 坐标区调整
set(gca, 'Box', 'off', ...                                                          % 边框
         'LineWidth', 1, 'GridLineStyle', '-',...                                   % 坐标轴线宽
         'XGrid', 'off', 'YGrid', 'off','ZGrid', 'on', ...                          % 网格
         'TickDir', 'out', 'TickLength', [.015 .015], ...                           % 刻度
         'XMinorTick', 'off', 'YMinorTick', 'off',  'ZMinorTick', 'off',...         % 小刻度
         'XColor', [.1 .1 .1],  'YColor', [.1 .1 .1], 'ZColor', [.1 .1 .1])      % 坐标轴颜色
                                

% 字体和字号
set(gca, 'FontName', 'Times New Roman')
set([hXLabel, hYLabel,hZLabel], 'FontName', '宋体')
set(gca, 'FontSize', 10)
set([hXLabel, hYLabel,hZLabel], 'FontSize', 10)
set(hTitle, 'FontName', '宋体','FontSize', 10, 'FontWeight' , 'bold')

% 绘制各储能输出功率
figure(2)
bar3(cell2mat(Pes_g_total),0.6);
colormap(slanCM(100))
hXLabel = xlabel('节点');
hYLabel = ylabel('场景');
hZLabel = zlabel('输出功率p.u.');

hTitle=title('不同扰动场景下网侧储能输出功率');

% 坐标区调整
set(gca, 'Box', 'off', ...                                                          % 边框
         'LineWidth', 1, 'GridLineStyle', '-',...                                   % 坐标轴线宽
         'XGrid', 'off', 'YGrid', 'off','ZGrid', 'on', ...                          % 网格
         'TickDir', 'out', 'TickLength', [.015 .015], ...                           % 刻度
         'XMinorTick', 'off', 'YMinorTick', 'off',  'ZMinorTick', 'off',...         % 小刻度
         'XColor', [.1 .1 .1],  'YColor', [.1 .1 .1], 'ZColor', [.1 .1 .1],...      % 坐标轴颜色
         'Xticklabel',NodeESCandidate)      
                                

% 字体和字号
set(gca, 'FontName', 'Times New Roman')
set([hXLabel, hYLabel,hZLabel], 'FontName', '宋体')
set(gca, 'FontSize', 10)
set([hXLabel, hYLabel,hZLabel], 'FontSize', 10)
set(hTitle, 'FontName', '宋体','FontSize', 10, 'FontWeight' , 'bold')

figure(3)
bar3(cell2mat(Psg_total),0.6);
colormap(slanCM(100))
hXLabel = xlabel('节点');
hYLabel = ylabel('场景');
hZLabel = zlabel('输出功率p.u.');

hTitle=title('不同扰动场景下发电机输出功率');

% 坐标区调整
set(gca, 'Box', 'off', ...                                                          % 边框
         'LineWidth', 1, 'GridLineStyle', '-',...                                   % 坐标轴线宽
         'XGrid', 'off', 'YGrid', 'off','ZGrid', 'on', ...                          % 网格
         'TickDir', 'out', 'TickLength', [.015 .015], ...                           % 刻度
         'XMinorTick', 'off', 'YMinorTick', 'off',  'ZMinorTick', 'off',...         % 小刻度
         'XColor', [.1 .1 .1],  'YColor', [.1 .1 .1], 'ZColor', [.1 .1 .1])      % 坐标轴颜色
                                

% 字体和字号
set(gca, 'FontName', 'Times New Roman')
set([hXLabel, hYLabel,hZLabel], 'FontName', '宋体')
set(gca, 'FontSize', 10)
set([hXLabel, hYLabel,hZLabel], 'FontSize', 10)
set(hTitle, 'FontName', '宋体','FontSize', 10, 'FontWeight' , 'bold')


figure(4)
% 储能配置惯性
subplot(2,1,1)
b=bar(Jes_g.*Ses_g);
b.FaceColor='flat';
colortemp=slanCM(188);
b.CData = colortemp(255);
% 修饰一下
hXLabel = xlabel('节点');
hYLabel = ylabel('惯性(s)');
set(gca, 'FontName', 'Times New Roman')
set([hXLabel, hYLabel], 'FontName', '宋体')
ax=gca;hold on;grid on
ax.LineWidth=1.2;
ax.XMinorTick='on';
ax.YMinorTick='on';
ax.ZMinorTick='on';
ax.GridLineStyle=':';
% 储能配置容量
subplot(2,1,2)
b=bar(Ses_g);
b.FaceColor='flat';
colortemp=slanCM(188);
b.CData = colortemp(255);
% 修饰一下
hXLabel = xlabel('节点');
hYLabel = ylabel('容量(p.u.)');
set(gca, 'FontName', 'Times New Roman')
set([hXLabel, hYLabel], 'FontName', '宋体')
ax=gca;hold on;grid on
ax.LineWidth=1.2;
ax.XMinorTick='on';
ax.YMinorTick='on';
ax.ZMinorTick='on';
ax.GridLineStyle=':';