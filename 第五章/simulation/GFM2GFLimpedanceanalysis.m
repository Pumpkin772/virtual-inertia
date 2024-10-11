clear
clc
close all
% 参考非对角元素影响word文档
W0 = 314.1593;
LF = 0.05;% pu
CF = 0.06;% 0.06pu
Kccp = 0.3;
Kcci = 10;
Kvf = 1;
Tvf = 0.02;% 0.02
Kvcp = 2;
Kvci = 10;
J = 0;
D = 100;


Kpllp = 5;
Kplli = 9500;

%% 
syms s
% GFM参数
% dq阻抗 不同Tvf
figure(1)
op = bodeoptions;
op.FreqUnits = 'Hz';
op.PhaseWrapping = 'on';
op.MagUnits = 'abs';
op.FreqScale = 'linear';
op.Title.FontSize = 10;
op.XLabel.FontSize = 10;
op.YLabel.FontSize = 10;
op.Xlim=[1 100];
PIvc = Kvcp + Kvci/s;
PIcc = Kccp + Kcci/s;
% Z0_S_simp = LF/W0/(PIcc*(PIvc+s*LF/W0)); % 这里除了S
% [num,den]=numden(Z0_S_simp);
% Num=sym2poly(num);  % 返回多项式项式系数
% Den=sym2poly(den);
% tf_Z0_S_simp = tf(Num,Den);
% bode(tf_Z0_S_simp*W0,op);hold on

Tvfk = [0.02 0.01 0.005 0];
for k = 1:length(Tvfk)
fvf = Kvf/(Tvfk(k)*s+1);
GI = PIcc/(PIcc+s*LF/W0);
YVF = (1-fvf)/(PIcc + s*LF/W0);
Y0 = (YVF + GI*PIvc+GI*s*CF/W0)/(1-GI);
Y0 = (1-fvf+PIcc*(PIvc+s*CF/W0))/(s*LF/W0);
Z0 = 1/Y0;
Z0_S = Z0/s;

[num,den]=numden(Z0_S);
Num=sym2poly(num);  % 返回多项式项式系数
Den=sym2poly(den);
tf_Z0_S = tf(Num,Den);

bode(tf_Z0_S*W0,op);hold on
end
h = findobj(gcf,'type','line');
set(h,'LineWidth',1)
h = findobj(gcf,'type','Axes');
set(h,'xtick',0:20:100);
set(h,'GridLineStyle',':','GridColor','[0.5 0.5 0.5]','GridAlpha',1);
set(h,'LineWidth',1);
set(h,'FontName','Times New Roman','FontSize',10);
grid on
% legend要在最下面
set(gcf,'position',[50 50 350 300]);
hleg = legend('{\it{T}}_{vf}=0.02','{\it{T}}_{vf}=0.01','{\it{T}}_{vf}=0.005','{\it{T}}_{vf}=0');
hleg = legend('{\it{T}}_{vf}=0.02','{\it{T}}_{vf}=0.01','{\it{T}}_{vf}=0.005','{\it{T}}_{vf}=0');% 这里要执行两次
set(hleg,'FontName','Times New Roman','FontSize',8);
hleg.NumColumns = 2;

%% 

op.Xlim=[10 100];
op.Title.String = 'J=0';
op.MagUnits = 'dB';
Vgfm0 = 1;
Id0 = 1;
theta0GFM = pi/4;


J=0;
Dk=[1 10 50 100];
for k=1:length(Dk)
Zqq=Vgfm0*(Id0/Y0+Vgfm0)/((J*s^2+Dk(k)*s)/W0)*sin(theta0GFM)*cos(theta0GFM);
figure(2)
Zqq_S = Zqq/s;
[num,den]=numden(Zqq_S);
Num=sym2poly(num);  % 返回多项式项式系数
Den=sym2poly(den);
tf_Zqq_S = tf(Num,Den);
bode(tf_Zqq_S*W0,op);hold on
end
h = findobj(gcf,'type','line');
set(h,'LineWidth',1)
h = findobj(gcf,'type','Axes');
set(h,'xtick',0:20:100);
set(h,'GridLineStyle',':','GridColor','[0.5 0.5 0.5]','GridAlpha',1);
set(h,'LineWidth',1);
set(h,'FontName','Times New Roman','FontSize',10);
grid on
% legend要在最下面
set(gcf,'position',[50 50 350 300]);
hleg = legend('{\it{D}}=1','{\it{D}}=10','{\it{D}}=50','{\it{D}}=100');
hleg = legend('{\it{D}}=1','{\it{D}}=10','{\it{D}}=50','{\it{D}}=100');% 这里要执行两次
set(hleg,'FontName','Times New Roman','FontSize',8);
hleg.NumColumns = 2;
%% 
op.Title.String = 'D=0';
D=0;
Jk=[0.1 1 5 10];
for k=1:length(Jk)
Zqq=Vgfm0*(Id0/Y0+Vgfm0)/((Jk(k)*s^2+D*s)/W0)*sin(theta0GFM)*cos(theta0GFM);
figure(3)
Zqq_S = Zqq/s;
[num,den]=numden(Zqq_S);
Num=sym2poly(num);  % 返回多项式项式系数
Den=sym2poly(den);
tf_Zqq_S = tf(Num,Den);
bode(tf_Zqq_S*W0,op);hold on
end
h = findobj(gcf,'type','line');
set(h,'LineWidth',1)
h = findobj(gcf,'type','Axes');
set(h,'xtick',0:20:100);
set(h,'GridLineStyle',':','GridColor','[0.5 0.5 0.5]','GridAlpha',1);
set(h,'LineWidth',1);
set(h,'FontName','Times New Roman','FontSize',10);
grid on
% legend要在最下面
set(gcf,'position',[50 50 350 300]);
hleg = legend('{\it{J}}=0.1','{\it{J}}=1','{\it{J}}=5','{\it{J}}=10');
hleg = legend('{\it{J}}=0.1','{\it{J}}=1','{\it{J}}=5','{\it{J}}=10');% 这里要执行两次
set(hleg,'FontName','Times New Roman','FontSize',8);
hleg.NumColumns = 2;




