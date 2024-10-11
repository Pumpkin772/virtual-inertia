clear
clc

RGB = [217 178 172;161 63 48; 
       251 206 189; 246 133 90;
       176 219 176;40 158 40; 
       190 174 211; 147 102 189;
       204 228 241;0 119 187;
       255 240 174;255 217 48;];%
   
figure(1)
configureresult = zeros(5,1);
load('twoGFMconfig');
configureresult(1)=Objective;
load('fourGFMconfig');
configureresult(2)=Objective;
load('sixGFMconfig');
configureresult(3)=Objective;
load('eightGFMconfig');
configureresult(4)=Objective;
load('tenGFMconfig');
configureresult(5)=Objective;

set(gcf,'unit','centimeters','position',[10 10 9 5]);
plot(2:2:10,configureresult,'-s','MarkerSize',10,'MarkerFaceColor',RGB(2,:)./255,'MarkerEdgeColor','none','Color',RGB(2,:)./255,'LineWidth',1.5);hold on
set(gca,'linewidth',1,'fontsize',10,'fontname','Times');
grid on
set(gca,'GridLineStyle',':','GridColor','[0.5 0.5 0.5]','GridAlpha',1);   
xlabel('储能个数');
ylabel('成本');
