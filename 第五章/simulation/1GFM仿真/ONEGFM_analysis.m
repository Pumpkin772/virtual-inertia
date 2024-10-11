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
s1=tf('s');

ws=2*pi*[1:1:100];
nu=length(ws);
s3=freqresp(s1,ws);

% GFM参数
% dq阻抗 不同Tvf
figure(1)
% op = bodeoptions;
% op.FreqUnits = 'Hz';
% op.PhaseWrapping = 'on';
% op.MagUnits = 'abs';
% op.FreqScale = 'linear';
% op.Title.FontSize = 10;
% op.XLabel.FontSize = 10;
% op.YLabel.FontSize = 10;
% op.Xlim=[1 100];
% PIvc = Kvcp + Kvci/s;
% PIcc = Kccp + Kcci/s;
% % Z0_S_simp = LF/W0/(PIcc*(PIvc+s*LF/W0)); % 这里除了S
% % [num,den]=numden(Z0_S_simp);
% % Num=sym2poly(num);  % 返回多项式项式系数
% % Den=sym2poly(den);
% % tf_Z0_S_simp = tf(Num,Den);
% % bode(tf_Z0_S_simp*W0,op);hold on

Tvfk = [0.02 0.01 0.005 0];
for k = 1:length(Tvfk)
    for ii=1:nu
        s = s3(ii);
        PIvc = Kvcp + Kvci/s;
PIcc = Kccp + Kcci/s;
fvf = Kvf/(Tvfk(k)*s+1);
GI = PIcc/(PIcc+s*LF/W0);
YVF = (1-fvf)/(PIcc + s*LF/W0);
Y0 = (YVF + GI*PIvc+GI*s*CF/W0)/(1-GI);
Y0 = (1-fvf+PIcc*(PIvc+s*CF/W0))/(s*LF/W0);
Z0 = 1/Y0;
% Z0_S = Z0/s;

% [num,den]=numden(Z0_S);
% Num=sym2poly(num);  % 返回多项式项式系数
% Den=sym2poly(den);
% tf_Z0_S = tf(Num,Den);

% bode(tf_Z0_S*W0,op);hold on
    R0(ii)=real(Z0);
    end
    plot(ws/2/pi,R0);hold on
end
ylim([-0.01 0.02]);
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