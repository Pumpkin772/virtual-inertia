clear
clc

RGB = [217 178 172;161 63 48; 
       251 206 189; 246 133 90;
       176 219 176;40 158 40; 
       190 174 211; 147 102 189;
       204 228 241;0 119 187];%255 240 174;255 217 48;
   
%% 数据处理
RoCoF_sim = zeros(8,10);
rang = 1:1000;
load('仿真结果/GFM/load1.mat')

RoCoF_sim(1,:)=min(out.RoCoF0.signals.values);
load('仿真结果/GFM/load2.mat')
rang = length(out.tout);
rang = 1:rang;
RoCoF_sim(2,:)=min(out.RoCoF0.signals.values);
load('仿真结果/GFM/load3.mat')
rang = length(out.tout);
rang = 1:rang;
RoCoF_sim(3,:)=min(out.RoCoF0.signals.values);
load('仿真结果/GFM/load4.mat')
rang = length(out.tout);
rang = 1:rang;
RoCoF_sim(4,:)=min(out.RoCoF0.signals.values);
load('仿真结果/GFM/load5.mat')
RoCoF_sim(5,:)=min(out.RoCoF0.signals.values);
load('仿真结果/GFM/load6.mat')
RoCoF_sim(6,:)=min(out.RoCoF0.signals.values);
load('仿真结果/GFM/load7.mat')
RoCoF_sim(7,:)=min(out.RoCoF0.signals.values);
load('仿真结果/GFM/load8.mat')
RoCoF_sim(8,:)=min(out.RoCoF0.signals.values);
RoCoF_sim(8,9)=-0.5235;
% load('lilunzhi.mat')
% 

x = 1:10;
y = mean(RoCoF_sim);

figure(1)
load('twoGFMpara.mat');
set(gcf,'unit','centimeters','position',[10 10 9 5]);
GO=bar(Jes_g.*Ses_g,1,'EdgeColor','k');hold off
GO(1).FaceColor = RGB(2,:)./255;
% lgd = legend({'仿真值','计算值'},'FontName','宋体','FontSize',10,'Location','southeast','NumColumns',2);
% errorbar(y,std(RoCoF-RoCoF_sim),'k','Linestyle','None','LineWidth', 1.2);
set(gca,'linewidth',1,'fontsize',10,'fontname','Times');
xlabel('节点','FontName','宋体','FontSize',10);
ylabel('等效惯性(s)','FontName','宋体','FontSize',10);
grid on;
set(gca,'GridLineStyle',':','GridColor','[0.5 0.5 0.5]','GridAlpha',1);

figure(2)
load('twoGFMpara.mat');
set(gcf,'unit','centimeters','position',[10 10 9 5]);
GO=bar(Ses_g*1000,1,'EdgeColor','k');hold off
GO(1).FaceColor = RGB(2,:)./255;
% lgd = legend({'仿真值','计算值'},'FontName','宋体','FontSize',10,'Location','southeast','NumColumns',2);
% errorbar(y,std(RoCoF-RoCoF_sim),'k','Linestyle','None','LineWidth', 1.2);
set(gca,'linewidth',1,'fontsize',10,'fontname','Times');
xlabel('节点','FontName','宋体','FontSize',10);
ylabel('功率容量(MW)','FontName','宋体','FontSize',10);
grid on;
set(gca,'GridLineStyle',':','GridColor','[0.5 0.5 0.5]','GridAlpha',1);

figure(3)
load('twoGFMconfig.mat');
set(gcf,'unit','centimeters','position',[10 10 9 5]);
GO=bar(Ses_g*1000,1,'EdgeColor','k');hold off
GO(1).FaceColor = RGB(2,:)./255;
% lgd = legend({'仿真值','计算值'},'FontName','宋体','FontSize',10,'Location','southeast','NumColumns',2);
% errorbar(y,std(RoCoF-RoCoF_sim),'k','Linestyle','None','LineWidth', 1.2);
set(gca,'linewidth',1,'fontsize',10,'fontname','Times');
xlabel('节点','FontName','宋体','FontSize',10);
ylabel('功率容量(MW)','FontName','宋体','FontSize',10);
grid on;
set(gca,'GridLineStyle',':','GridColor','[0.5 0.5 0.5]','GridAlpha',1);


figure(3)
set(gcf,'unit','centimeters','position',[10 10 9 5]);
load('仿真结果/GFM/load1.mat')
yanse=slanCM('tab10',10);
% yanse=slanCM(182,10);
p=plot(out.tout,out.Freq0.signals.values,'linewidth',1);
set(gca,'linewidth',1,'fontsize',10,'fontname','Times');
for i = 1:10
    p(i).Color=yanse(i,:);
end
xlabel('时间(s)','FontName','宋体','FontSize',10);
ylabel('频率(Hz)','FontName','宋体','FontSize',10);
grid on;
set(gca,'GridLineStyle',':','GridColor','[0.5 0.5 0.5]','GridAlpha',1);

figure(4)
set(gcf,'unit','centimeters','position',[10 10 9 5]);
load('仿真结果/GFM/load1.mat')
yanse=slanCM('tab10',10);
% yanse=slanCM(182,10);
p=plot(out.tout,out.RoCoF0.signals.values,'linewidth',1);
set(gca,'linewidth',1,'fontsize',10,'fontname','Times');
for i = 1:10
    p(i).Color=yanse(i,:);
end
xlabel('时间(s)','FontName','宋体','FontSize',10);
ylabel('频率变化率(Hz/s)','FontName','宋体','FontSize',10);
grid on;
set(gca,'GridLineStyle',':','GridColor','[0.5 0.5 0.5]','GridAlpha',1);


figure(5)
set(gcf,'unit','centimeters','position',[10 10 9 5]);
load('仿真结果/GFL/load1.mat')
yanse=slanCM('tab10',10);
% yanse=slanCM(182,10);
p=plot(out.tout,out.Freq0.signals.values,'linewidth',1);
set(gca,'linewidth',1,'fontsize',10,'fontname','Times');
for i = 1:10
    p(i).Color=yanse(i,:);
end
xlabel('时间(s)','FontName','宋体','FontSize',10);
ylabel('频率(Hz)','FontName','宋体','FontSize',10);
grid on;
set(gca,'GridLineStyle',':','GridColor','[0.5 0.5 0.5]','GridAlpha',1);

figure(6)
set(gcf,'unit','centimeters','position',[10 10 9 5]);
load('仿真结果/GFL/load1.mat')
yanse=slanCM('tab10',10);
% yanse=slanCM(182,10);
p=plot(out.tout,out.RoCoF0.signals.values,'linewidth',1);
set(gca,'linewidth',1,'fontsize',10,'fontname','Times');
for i = 1:10
    p(i).Color=yanse(i,:);
end
xlabel('时间(s)','FontName','宋体','FontSize',10);
ylabel('频率变化率(Hz/s)','FontName','宋体','FontSize',10);
grid on;
set(gca,'GridLineStyle',':','GridColor','[0.5 0.5 0.5]','GridAlpha',1);

figure(7)
set(gcf,'unit','centimeters','position',[10 10 9 5]);
GO=bar(RoCoF_sim',1,'EdgeColor','k');hold off
yanse=slanCM('Dark2',8);
set(gca,'linewidth',1,'fontsize',10,'fontname','Times');
for i = 1:8
    GO(i).FaceColor=yanse(i,:);
end
xlabel('时间(s)','FontName','宋体','FontSize',10);
ylabel('频率变化率(Hz/s)','FontName','宋体','FontSize',10);
grid on;
set(gca,'GridLineStyle',':','GridColor','[0.5 0.5 0.5]','GridAlpha',1);
% 
% figure(5)
% bar3(RoCoF_sim,0.6);

figure(8)
set(gcf,'unit','centimeters','position',[10 10 9 5]);
load('仿真结果/只计初始频率变化率/load1.mat')
yanse=slanCM('tab10',10);
% yanse=slanCM(182,10);
p=plot(out.tout,out.Freq0.signals.values,'linewidth',1);
set(gca,'linewidth',1,'fontsize',10,'fontname','Times');
for i = 1:10
    p(i).Color=yanse(i,:);
end
xlabel('时间(s)','FontName','宋体','FontSize',10);
ylabel('频率(Hz)','FontName','宋体','FontSize',10);
grid on;
set(gca,'GridLineStyle',':','GridColor','[0.5 0.5 0.5]','GridAlpha',1);

figure(9)
set(gcf,'unit','centimeters','position',[10 10 9 5]);
load('仿真结果/只计初始频率变化率/load1.mat')
yanse=slanCM('tab10',10);
% yanse=slanCM(182,10);
p=plot(out.tout,out.RoCoF0.signals.values,'linewidth',1);
set(gca,'linewidth',1,'fontsize',10,'fontname','Times');
for i = 1:10
    p(i).Color=yanse(i,:);
end
xlabel('时间(s)','FontName','宋体','FontSize',10);
ylabel('频率变化率(Hz/s)','FontName','宋体','FontSize',10);
grid on;
set(gca,'GridLineStyle',':','GridColor','[0.5 0.5 0.5]','GridAlpha',1);

%% 配置4个储能

figure(10)
set(gcf,'unit','centimeters','position',[10 10 9 5]);
load('仿真结果/四个GFM/load1.mat')
yanse=slanCM('tab10',10);
% yanse=slanCM(182,10);
p=plot(out.tout,out.RoCoF0.signals.values,'linewidth',1);
set(gca,'linewidth',1,'fontsize',10,'fontname','Times');
for i = 1:10
    p(i).Color=yanse(i,:);
end
xlabel('时间(s)','FontName','宋体','FontSize',10);
ylabel('频率变化率(Hz/s)','FontName','宋体','FontSize',10);
grid on;
set(gca,'GridLineStyle',':','GridColor','[0.5 0.5 0.5]','GridAlpha',1);

figure(11)
set(gcf,'unit','centimeters','position',[10 10 9 5]);
load('仿真结果/四个GFM/load1.mat')
yanse=slanCM('tab10',10);
% yanse=slanCM(182,10);
p=plot(out.tout,out.Freq0.signals.values,'linewidth',1);
set(gca,'linewidth',1,'fontsize',10,'fontname','Times');
for i = 1:10
    p(i).Color=yanse(i,:);
end
xlabel('时间(s)','FontName','宋体','FontSize',10);
ylabel('频率(Hz)','FontName','宋体','FontSize',10);
grid on;
set(gca,'GridLineStyle',':','GridColor','[0.5 0.5 0.5]','GridAlpha',1);

figure(12)
set(gcf,'unit','centimeters','position',[10 10 9 5]);
load('仿真结果/四个GFM/load1.mat')
yanse=slanCM('tab10',4);
% yanse=slanCM(182,10);

for i = 1:4
    p=plot(out.tout,out.PGFM.signals(i).values*1000,'linewidth',1);hold on
    set(gca,'linewidth',1,'fontsize',10,'fontname','Times');
    p.Color=yanse(i,:);
end
xlabel('时间(s)','FontName','宋体','FontSize',10);
ylabel('输出功率(Hz/s)','FontName','宋体','FontSize',10);
grid on;
set(gca,'GridLineStyle',':','GridColor','[0.5 0.5 0.5]','GridAlpha',1);