syms s
W0 = 2*pi*50;
Jii = 2*logspace(-2,1,200);
Dii = 2*logspace(-2,2,200);

theta0 = 0;
Vd0 = 1;
SGFM = 1;
LG = 0.06;
RG = 0;

XL = 0.3;
[JJ,DD]=meshgrid(Jii,Dii);
XL = XL + 0.06;
for SGFM=0.1:0.1:1

% stability=zeros(20,20);
% for i = 1 : 20
%     for j = 1:20
%         J = JJ(i,j);
%         D = DD(i,j);
%         Zmtheta = -Vd0*(Vd0)/((J*s^2+D*s)/W0)*[sin(theta0)*cos(theta0) -sin(theta0)^2;-cos(theta0)^2 -sin(theta0)*cos(theta0)];
%         Y0 = 1/(0.002)*SGFM;
%         Zxygfm = diag([inv(Y0) inv(Y0)])+Zmtheta+LG*[s/W0 -1;1 s/W0]+RG*eye(2);
%         Q = 1/XL;
%         tao = 0;
%         gama = 1/((s+tao)^2/W0+W0)*[s+tao W0;-W0 s+tao];
%         eigvec = double(solve(det(inv(Zxygfm)+kron(Q,gama)),s));
%         if eigvec<0
%             stability(i,j) = 1;
%         else
%             stability(i,j) = 0;
%         end
%     end
% end


tao = 0.002/SGFM/(XL);
stability2=((tao^2+W0^2)*DD.*((2*tao*JJ+DD).*(JJ*(tao^2+W0^2)+2*tao*DD)-JJ*(tao^2+W0^2).*DD)./(2*tao*JJ+DD)-(2*tao*JJ+DD)*W0^2*Vd0^2/XL)>0;

P = [JJ(stability2),DD(stability2)];
k = boundary(P);
% plot(JJ(stability2),DD(stability2),'*');
plot(P(k,1),P(k,2));hold on
end


