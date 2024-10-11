% 多GFL和多GFM的惯性、下垂以及有功功率优化
clear
clc
close all
%% 
% casename = 'fixedGFMJ_D';
casename = 'varGFMJ_D';
% casename = 'fixedGFMJ_D_low';
% casename = 'withoutsmallsignalconstraint';
% casename = 'test';
% casename = 'varGFMJ_D_fixed';
TIME = 24;
marginlist = [0.1;0.5];
for mar = 1 : 2
    margin = marginlist(mar);
% 新能源预测输出功率
% PGFL_forecast = [0.23,0.23,0.23,0.23,0.23,0.23,0.23,0.23,0.59,0.59,1,1,0.59,0.59,1,1,1,1,1,0.23,0.23,0.23,0.23,0.23];
% PGFL_forecast = kron(PGFL_forecast,ones(9,1));
winddata = load('nrel_wind_power.csv');
winddata = reshape(winddata(1,2:end),[96,60])';
windsum = sum(winddata,2);
[B,I]=sort(windsum,'descend');
PGFL_forecast = winddata(I([1:5,21:24]),:);

if strcmp(casename,'test')
    PGFL_forecast = winddata([1 4 5 7 8 10 11 13 14],:);
end
if TIME~=96
temp = zeros(9,24);
for tim = 1:24
    temp(:,tim) = mean(PGFL_forecast(:,4*(tim-1)+1:4*tim),2);
end
PGFL_forecast=temp;
end
% GFM的强度
load('deltaGFM_strength');
% GFM的J、D稳定范围
load('GFMSTABLE_J_D');
% 参数
% Kpllp_list = [8,8,8,8,8,8,8,8,8];
Kpllp_list = [8,8,8,8,8,8.5,8.5,8.5,8.5];
Kplli_list = [9500,9500,9300,9800,9600,9300,9200,9500,9500];
if strcmp(casename,'test')
    Kpllp_list = [8,8,8,8,8,7,7,7,7];
    Kplli_list = [9500,9500,9300,9800,9600,9300,9200,9500,9500];
end
W0 = 2*pi*50;
Kp_avg = Kpllp_list*pp_2;
Ki_avg = Kplli_list*pp_2;
tao = 0;
Uavg = 1;

% bound=W0*Kp_avg*Uavg*lamb2;% PLL小扰动约束
%% 
% 新能源数据
NGFL = size(PGFL_forecast,1); % GFL新能源数量
period = size(PGFL_forecast,2); % 时间段

% GFM储能数据
PGFM_max = [0.8;1.2]; % GFM储能最大输出功率约束
NGFM = length(PGFM_max); % GFM储能数量
EGFM_cap = 2*PGFM_max; % GFM储能的能量容量
EGFM_ini = 0.5*EGFM_cap; % GFM储能初始容量0.2
EGFM_min = 0.2*EGFM_cap; % GFM储能最小容量0.2
EGFM_max = 0.95*EGFM_cap; % GFM储能最大容量0.95
effi = 0.98; % 储能充放电效率

% 频率限制数据
RoCoF_max = 2;% 最大频率变化率Hz/s
freq_max  = 0.5;% 最大频率偏差Hz
freq_ss   = 0.25;% 稳态频率偏差Hz

% 惯性、下垂系数
JGFMmin = 0.01;% GFM不能为0
Jmax = 100;% 最大惯性p.u.
Dmax = 100;% 最大下垂p.u.
dis  = 0.1;% 设置扰动为5%的输出功率
% GFM的J、D稳定范围



%% 决策变量
PGFL = sdpvar(NGFL,period,'full'); % 新能源输出功率
DGFL = sdpvar(NGFL,period,'full'); % 新能源提供的阻尼D
ReGFL = sdpvar(NGFL,period,'full'); % 新能源备用功率容量
P_Dmax = sdpvar(NGFL,period,'full'); % 新能源小扰动最大输出功率
% Ytemp = sdpvar(NGFL,period,'full'); % 

PGFM_c = sdpvar(NGFM, period, 'full'); % 储能充电功率
PGFM_d = sdpvar(NGFM, period, 'full'); % 储能放电功率
IGFM = binvar(NGFM, period, 'full'); % 储能充放电状态，取1时表示充电，取0时表示放电
EGFM = sdpvar(NGFM, period, 'full'); % 储能存储电量
LGFM = sdpvar(NGFM, period, 'full'); % 储能运行损耗
JGFM = sdpvar(NGFM, period, 'full'); % 储能提供的惯性J
JGFM0 = sdpvar(1, period, 'full'); % 储能提供的惯性J0
DGFM0 = sdpvar(1, period, 'full'); % 储能提供的惯性D0
DGFM = sdpvar(NGFM, period, 'full'); % 储能提供的阻尼D
ReGFM = sdpvar(NGFM,period,'full'); % 储能备用功率容量
deltaYgGFM = sdpvar(1,period,'full'); % 储能输出阻抗

Jfixed = sdpvar(1, 1, 'full'); % 储能提供的惯性J0
Dfixed = sdpvar(1, 1, 'full'); % 储能提供的惯性D0

Pload = sdpvar(1,period,'full'); % 所提供的负荷功率
%% 约束条件
cons = [];

% 储能充放电最大充放电功率
for tim = 1:NGFM
cons = [cons,...
    0 <= PGFM_c(tim,:) <= PGFM_max(tim)*IGFM(tim,:)];
cons = [cons,...
    0 <= PGFM_d(tim,:) <= PGFM_max(tim)*(1-IGFM(tim,:))];
end
% 储能剩余电量更新
cons = [cons,...
    EGFM(:,1) == EGFM_ini + effi*PGFM_c(:,1) - PGFM_d(:,1)/effi];
for tim = 2:period
cons = [cons,...
    EGFM(:,tim) == EGFM(:,tim-1) + effi*PGFM_c(:,tim) - PGFM_d(:,tim)/effi];
end
% 储能剩余电量约束
cons = [cons,...
    EGFM_min.*ones(NGFM, period) <= EGFM,...
    EGFM <= EGFM_max.*ones(NGFM, period)];
cons = [cons,...
    EGFM(:,period) == EGFM_ini];
% 储能损耗
cons = [cons,...
    LGFM >= (1/effi - 1)*PGFM_d, ...
    LGFM >= (1 - effi)*PGFM_c];

% 新能源输出功率约束
cons = [cons,0 <= PGFL,...
    PGFL <= PGFL_forecast];

% 功率平衡约束
for tim = 1:period
cons = [cons,Pload(tim) == sum(PGFM_d(:,tim)) - sum(PGFM_c(:,tim)) + sum(PGFL(:,tim))];
end

% 频率约束（这里还没有考虑频率最低点）
cons = [cons,...
    Pload*dis*50/RoCoF_max <= ones(1,NGFM)*JGFM,...
    Pload*dis*50/freq_ss <= 0*ones(1,NGFL)*DGFL + ones(1,NGFM)*DGFM];

% 备用容量约束
cons = [cons,...
    DGFL*freq_ss/50 == ReGFL,... % 确保能够满足频率调节
    JGFM*RoCoF_max/50 + DGFM*freq_ss/50 == ReGFM,...
    ReGFL + PGFL <= PGFL_forecast,... % 不能超过预测新能源功率
    PGFM_d - PGFM_c + ReGFM <= PGFM_max.*ones(NGFM, period)
    ];

% 控制器参数约束
cons = [cons,...
    JGFMmin <= JGFM0,...
    JGFM0 <= Jmax,...
    0 <= DGFL,...
    DGFL<= Dmax,...
    0 <= DGFM0,...
    DGFM0<= Dmax];% 12.5恒定惯性、下垂系数
if strcmp(casename,'fixedGFMJ_D')
    cons = [cons,...
    JGFM0==12.5*ones(1,period),...
    DGFM0==100*ones(1,period)];% 12.5恒定惯性、下垂系数
end
if strcmp(casename,'fixedGFMJ_D_low')
    cons = [cons,...
    JGFM0==5*ones(1,period),...
    DGFM0==50*ones(1,period)];% 5恒定惯性、50下垂系数
end
if strcmp(casename,'varGFMJ_D_fixed')
    cons = [cons,...
    JGFM0==Jfixed*ones(1,period),...
    DGFM0==Dfixed*ones(1,period)];% 5恒定惯性、50下垂系数
end
for tim = 1:NGFM
    cons = [cons,...
        JGFM(tim,:)==PGFM_max(tim)*JGFM0,...
        DGFM(tim,:)==PGFM_max(tim)*DGFM0];
end

% GFM的J、D范围约束

% 小扰动稳定约束（这里还没有添加）
cons = [cons,...
%     deltaYgGFM == interp2(JJ,DD,deltaGFM_strength,JGFM0,DGFM0,'sos2'),...% GFM输出阻抗与J的关系sos2，这里需要拟合的是Yg，因为并联
%     P_Dmax == interp2(YY,JJ,Pmax,YgGFM,JGFL,'sos2'),...% GFL输出功率与阻抗、J、D的关系，这里是YgGFM+1，假定原有的连接阻抗为1（还没有设置功率限制）
%     Ytemp == 2,...
%     PGFL - 0*DGFL <= P_Dmax,...
    ];
    for tim = 1:period
        cons = [cons,...
            %     Ki_avg*pp_2'*PGFL(:,i)<= (lamb2)*W0*Kp_avg*Uavg;
            deltaYgGFM(:,tim) == interp2(JJ,DD,deltaGFM_strength,JGFM0(:,tim),DGFM0(:,tim),'sos2');
            ];
    end
if ~strcmp(casename,'withoutsmallsignalconstraint')
    
    for tim = 1:period
        cons = [cons,...
            %     Ki_avg*pp_2'*PGFL(:,i)<= (lamb2)*W0*Kp_avg*Uavg;
            Ki_avg*pp_2'*PGFL(:,tim) + tao*Kp_avg*pp_2'*PGFL(:,tim)<= (lamb2-margin+deltaYgGFM(:,tim))*W0*Kp_avg*Uavg;
            ];
    end
end
%% 目标函数
obj = 0;
for tim = 1:period
   obj = obj +...
       ones(1,NGFL)*(PGFL_forecast(:,tim) - PGFL(:,tim)) ...% 弃风惩罚
       + 0.0*ones(1,NGFL)*ReGFL(:,tim)... % 新能源备用容量
       + 0.0*ones(1,NGFM)*ReGFM(:,tim) ... % 储能备用容量
       + 0.01*ones(1,NGFM)*LGFM(:,tim); % 储能损耗
end
%% 求解器设置
ops = sdpsettings;
ops.solver = 'gurobi';
ops.gurobi.TimeLimit = 500;
% ops.showprogress = 0;
% ops.verbose = 0;
% ops.cachesolvers = 1;
% ops.savesolveroutput = 1;
Result = optimize(cons, obj, ops);

% if Result.problem == 0
    disp(Result.info);
    Result.obj = value(obj);
    Result.PGFL = value(PGFL);
    Result.DGFL = value(DGFL);
    Result.ReGFL = value(ReGFL);
    Result.P_Dmax = value(P_Dmax);

    
    Result.PGFM_c = value(PGFM_c);
    Result.PGFM_d = value(PGFM_d);
    Result.EGFM = value(EGFM);
    Result.LGFM = value(LGFM);
    Result.JGFM = value(JGFM);
    Result.DGFM = value(DGFM);
    Result.ReGFM = value(ReGFM);
    Result.deltaYgGFM = value(deltaYgGFM);
    
    Result.P_Dmax = value(P_Dmax);
    
    Result.Pload = value(Pload);
    Result.PGFL_forecast = PGFL_forecast;
    Result.REScurtail = sum(PGFL_forecast,'all') - sum(Result.PGFL,'all');
    
    

% else
    disp(Result.info);
% end
%%
check = zeros(period,1);
for tim = 1:period
    check(tim) = Ki_avg*pp_2'*Result.PGFL(:,tim) + tao*Kp_avg*pp_2'*Result.PGFL(:,tim) - (lamb2+Result.deltaYgGFM(:,tim))*W0*Kp_avg*Uavg;
%     check(i) = Ki_avg*pp_2'*Result.PGFL(:,i) + tao*Kp_avg*pp_2'*Result.PGFL(:,i) - (lamb2)*W0*Kp_avg*Uavg;
end
%%
PGFM_max = [0.8;1.2]; % GFM储能最大输出功率约束
NGFM = length(PGFM_max); % GFM储能数量
EGFM_cap = 2*PGFM_max; % GFM储能的能量容量
dis  = 0.1;
NGFL = 9;
    figure(1)
%     stairs(PGFL_forecast,'LineWidth',1.5);hold on
subplot(2,1,1)
stairs(Result.PGFL_forecast','LineWidth',1.5);hold on
    stairs((Result.PGFL)','LineWidth',1.5);
    subplot(2,1,2)
    stairs(sum(Result.PGFL_forecast)','LineWidth',1.5);hold on
    stairs(Result.Pload,'LineWidth',1.5);
    xlabel('Time (h)');
ylabel('P (p.u.)');
set(gca,'GridLineStyle',':','GridColor','[0.5 0.5 0.5]','GridAlpha',1);
set(gca,'LineWidth',1);
set(gca,'FontName','Times New Roman','FontSize',10);
grid on
set(gcf,'position',[50 50 350 300]);
figure(2)
stairs((Result.PGFM_d-Result.PGFM_c)','LineWidth',1.5);
xlabel('Time (h)');
ylabel('P (p.u.)');
set(gca,'GridLineStyle',':','GridColor','[0.5 0.5 0.5]','GridAlpha',1);
set(gca,'LineWidth',1);
set(gca,'FontName','Times New Roman','FontSize',10);
grid on
set(gcf,'position',[50 50 350 300]);
figure(3)
stairs((Result.EGFM./EGFM_cap)','LineWidth',1.5);
xlabel('Time (h)');
ylabel('SoC');
set(gca,'GridLineStyle',':','GridColor','[0.5 0.5 0.5]','GridAlpha',1);
set(gca,'LineWidth',1);
set(gca,'FontName','Times New Roman','FontSize',10);
grid on
set(gcf,'position',[50 50 350 300]);
figure(4)
stairs((Result.JGFM)','LineWidth',1.5);
xlabel('Time (h)');
ylabel('Inertia (s)');
set(gca,'GridLineStyle',':','GridColor','[0.5 0.5 0.5]','GridAlpha',1);
set(gca,'LineWidth',1);
set(gca,'FontName','Times New Roman','FontSize',10);
grid on
set(gcf,'position',[50 50 350 300]);
figure(5)
stairs((Result.deltaYgGFM)','LineWidth',1.5);
xlabel('Time (h)');
ylabel('Strength (p.u.)');
set(gca,'GridLineStyle',':','GridColor','[0.5 0.5 0.5]','GridAlpha',1);
set(gca,'LineWidth',1);
set(gca,'FontName','Times New Roman','FontSize',10);
grid on
set(gcf,'position',[50 50 350 300]);
figure(6)
stairs((Result.DGFM)','LineWidth',1.5);
xlabel('Time (h)');
ylabel('Damping (p.u.)');
set(gca,'GridLineStyle',':','GridColor','[0.5 0.5 0.5]','GridAlpha',1);
set(gca,'LineWidth',1);
set(gca,'FontName','Times New Roman','FontSize',10);
grid on
set(gcf,'position',[50 50 350 300]);
figure(7)
stairs((Result.ReGFM)','LineWidth',1.5);
xlabel('Time (h)');
ylabel('Reserve (p.u.)');
set(gca,'GridLineStyle',':','GridColor','[0.5 0.5 0.5]','GridAlpha',1);
set(gca,'LineWidth',1);
set(gca,'FontName','Times New Roman','FontSize',10);
grid on
set(gcf,'position',[50 50 350 300]);
figure(8)
stairs((Result.Pload*dis*50./(ones(1,NGFM)*Result.JGFM))','LineWidth',1.5);
xlabel('Time (h)');
ylabel('RoCoF (Hz/s)');
set(gca,'GridLineStyle',':','GridColor','[0.5 0.5 0.5]','GridAlpha',1);
set(gca,'LineWidth',1);
set(gca,'FontName','Times New Roman','FontSize',10);
grid on
set(gcf,'position',[50 50 350 300]);
figure(9)
stairs((Result.Pload*dis*50./(ones(1,NGFL)*Result.DGFL + ones(1,NGFM)*Result.DGFM))','LineWidth',1.5);
xlabel('Time (h)');
ylabel('Steady State (Hz)');
set(gca,'GridLineStyle',':','GridColor','[0.5 0.5 0.5]','GridAlpha',1);
set(gca,'LineWidth',1);
set(gca,'FontName','Times New Roman','FontSize',10);
grid on
set(gcf,'position',[50 50 350 300]);
%% 
calculate = 1;
save(['OptimizeResult/',casename,'_',num2str(period),'_',num2str(margin*100)],'Result');
if calculate == 1
GFL_con = [39,2,3,5,7,8,9,10,11];
GFM_con = [4,6];
eigvalue = cell(period,1);
eigvalue_GFL = zeros(period,1);
eigvalue_GFL_simple_0 = zeros(period,1);
eigvalue_GFL_simple_1 = zeros(period,1);
for tim = 1:period
    file = xlsread('ModifiedIEEE39_11gen.xlsx','Bus');
    file(GFL_con,5) = Result.PGFL(:,tim);
    file(GFM_con,5) = Result.PGFM_d(:,tim)-Result.PGFM_c(:,tim);
    xlswrite('ModifiedIEEE39_11gen.xlsx',file,'Bus','A4:L42');
    file = xlsread('ModifiedIEEE39_11gen.xlsx','Apparatus');
    file = file(2:end,:);
    file(GFM_con,9) = Result.DGFM(:,tim)./PGFM_max;
    file(GFM_con,13) = Result.JGFM(:,tim)./PGFM_max;
    file([1,2,3,5,7,8,9,10,11],8) = Kpllp_list;
    file([1,2,3,5,7,8,9,10,11],10) = Kplli_list;
    xlswrite('ModifiedIEEE39_11gen.xlsx',file,'Apparatus','A5:N15');
    run('ModifiedIEEE39_11gen')
    eigvalue(tim) = {EigVec};
    EigVec = sort(EigVec(imag(EigVec)<110&imag(EigVec)>90,1),'descend','ComparisonMethod','real');
    eigvalue_GFL(tim) = EigVec(1);
    
    U = zeros(NGFL,1);
Theta = zeros(NGFL,1);
for k = 1:NGFL
    U(k) = PowerFlow{1, GFL_con(k)}(3);
    Theta(k) = PowerFlow{1, GFL_con(k)}(4);
end
Seqp = 0;
Seqi = 0;
Ueqp = 0;
Ueqi = 0;
for k = 1 : NGFL
    kpllp=Kpllp_list(k);
    kplli=Kplli_list(k);
    Ueqp = Ueqp + pp_2(k)*kpllp*U(k);
    Ueqi = Ueqi + pp_2(k)*kplli*U(k);
    for j = 1:NGFL
        Seqp = Seqp + kpllp*pp_2(k)*pp_2(j)/U(j)*(ListBus(GFL_con(j),5)*cos(Theta(k)-Theta(j))+1i*ListBus(GFL_con(j),5)*sin(Theta(k)-Theta(j)));
        Seqi = Seqi + kplli*pp_2(k)*pp_2(j)/U(j)*(ListBus(GFL_con(j),5)*cos(Theta(k)-Theta(j))+1i*ListBus(GFL_con(j),5)*sin(Theta(k)-Theta(j)));
    end
end
Peqp = real(Seqp);
Qeqp = imag(Seqp);
Peqi = real(Seqi);
Qeqi = imag(Seqi);
Aeq = zeros(2,2);
sigma1 = 1/(lamb2+Result.deltaYgGFM(:,tim));
wref = W0;
Aeq(1,1)=sigma1*Peqi/(wref-sigma1*Peqp);
Aeq(1,2)=((sigma1*tao*Peqi-wref*Ueqi+wref*sigma1*Qeqi)*(wref-sigma1*Peqp)+sigma1*Peqi*(-wref*Ueqp+wref*sigma1*Qeqp+sigma1*tao*Peqp))/(wref^2-wref*sigma1*Peqp);
Aeq(2,1)=wref/(wref-sigma1*Peqp);
Aeq(2,2)=(-wref*Ueqp+sigma1*wref*Qeqp+sigma1*tao*Peqp)/(wref-sigma1*Peqp);
lamb_Asimp=eig(Aeq);
eigvalue_GFL_simple_0(tim) = lamb_Asimp(imag(lamb_Asimp)>0);

% 近似后
Peqp = 0;
Qeqp = 0;
Peqi = 0;
Qeqi = 0;
kpllpeq = Kpllp_list*pp_2;
kpllieq = Kplli_list*pp_2;
Uavg = mean(U);
Ueqi = kpllieq*Uavg;
Ueqp = kpllpeq*Uavg;
for j = 1:NGFL
    Peqp = Peqp + kpllpeq*ListBus(GFL_con(j),5)*pp_2(j)/Uavg;
    Peqi = Peqi + kpllieq*ListBus(GFL_con(j),5)*pp_2(j)/Uavg;
end
Aeq(1,1)=sigma1*Peqi/(wref-sigma1*Peqp);
Aeq(1,2)=((sigma1*tao*Peqi-wref*Ueqi+wref*sigma1*Qeqi)*(wref-sigma1*Peqp)+sigma1*Peqi*(-wref*Ueqp+wref*sigma1*Qeqp+sigma1*tao*Peqp))/(wref^2-wref*sigma1*Peqp);
Aeq(2,1)=wref/(wref-sigma1*Peqp);
Aeq(2,2)=(-wref*Ueqp+sigma1*wref*Qeqp+sigma1*tao*Peqp)/(wref-sigma1*Peqp);
lamb_Asimp2=eig(Aeq);
eigvalue_GFL_simple_1(tim) = lamb_Asimp2(imag(lamb_Asimp2)>0);
end
save(['OptimizeResult/',casename,'_',num2str(period)],'Result','eigvalue','eigvalue_GFL','eigvalue_GFL_simple_0','eigvalue_GFL_simple_1');
% 简化模型
eigvalue_GFL_simple = zeros(period,1);
for tim = 1:period
Peqp = 0;
Qeqp = 0;
Peqi = 0;
Qeqi = 0;
kpllpeq = Kp_avg;
kpllieq = Kplli_list*pp_2;
Uavg = 1;
Ueqi = kpllieq*Uavg;
Ueqp = kpllpeq*Uavg;
for j = 1:NGFL
    Peqp = Peqp + kpllpeq*Result.PGFL(j,tim)*pp_2(j)/Uavg;
    Peqi = Peqi + kpllieq*Result.PGFL(j,tim)*pp_2(j)/Uavg;
end
sigma1 = 1/(lamb2+Result.deltaYgGFM(:,tim));
wref = W0;
Aeq(1,1)=sigma1*Peqi/(wref-sigma1*Peqp);
Aeq(1,2)=((sigma1*tao*Peqi-wref*Ueqi+wref*sigma1*Qeqi)*(wref-sigma1*Peqp)+sigma1*Peqi*(-wref*Ueqp+wref*sigma1*Qeqp+sigma1*tao*Peqp))/(wref^2-wref*sigma1*Peqp);
Aeq(2,1)=wref/(wref-sigma1*Peqp);
Aeq(2,2)=(-wref*Ueqp+sigma1*wref*Qeqp+sigma1*tao*Peqp)/(wref-sigma1*Peqp);
lamb_Asimp2=eig(Aeq);
eigvalue_GFL_simple(tim) = lamb_Asimp2(imag(lamb_Asimp2)>0);
end
save(['OptimizeResult/',casename,'_',num2str(period),'_',num2str(margin*100)],'Result','eigvalue','eigvalue_GFL','eigvalue_GFL_simple_0','eigvalue_GFL_simple_1','eigvalue_GFL_simple');
end
end