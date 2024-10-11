function [I_new,Y_new]=kron_user(Y,CZno,CIno,I,YL)
% |I1|   |Y11 Y12 ... Y1N| |V1|
% |I2| = |Y21 Y21 ... Y2N| |V2|
% |..|   |               | |V3|
% |IN| = |YN1 YN2 ... YNN| |VN|
% 消去恒阻抗、恒电流节点
% 恒阻抗节点：Ix = Yx1V1+Yx2V2+...+YxNVN
%            Ix = -YLxVx(电流注入方向)
%            0  =  Yx1V1+Yx2V2+...+(Yxx+YLx)Vx+...+YxNVN
%            Vx =  -Yx1/(Yxx+YLx)V1-...-YxN/(Yxx+YLx)VN
%            Ii = (Yi1-YixYx1/(Yxx+YLx))V1+...+(YiN-YixYxN/(Yxx+YLx))VN
% 恒电流节点：Ix = Yx1V1+Yx2V2+...+YxNVN
%             Vx =  -Yx1/(Yxx)*V1-...-YxN/(Yxx)*VN + 1/Yxx*Ix
%             这里增加了一个常数给每个节点
%             Ii = (Yi1-YixYx1/(Yxx))V1+...+(YiN-YixYxN/(Yxx))VN + Yix/Yxx*Ix
n0  = length(I);
nCZ = length(CZno);
nCI = length(CIno);
newno = setdiff(1:n0,[CZno,CIno]);
nnew = n0 - nCZ - nCI;
Y_new = zeros(nnew,nnew);
I_new = zeros(nnew,1);

reduction_mat = [];
Y1 = zeros(n0,n0);
I1 = zeros(n0,1);
% 消去恒电流节点，附加电流
for k = 1:nCI
    for i = 1:n0
        if(ismember(i,reduction_mat))
            continue;
        end
        for j = 1:n0
            if(ismember(j,reduction_mat))
                continue;
            end
            Y1(i,j) = Y(i,j) - Y(i,CIno(k))*Y(CIno(k),j)/Y(CIno(k),CIno(k));
            I1(i) = I(i) - Y(i,CIno(k))/Y(CIno(k),CIno(k))*I(CIno(k));
        end
    end
    reduction_mat = [reduction_mat CIno(k)];
    Y = Y1;
    I = I1;
end
Y_temp = Y;
% 消去恒阻抗节点以及恒电流节点对导纳矩阵影响
for k = 1:nCZ
    for i = 1:n0
        if(ismember(i,reduction_mat))
            continue;
        end
        for j = 1:n0
            if(ismember(j,reduction_mat))
                continue;
            end
            Y1(i,j) = Y(i,j) - Y(i,CZno(k))*Y(CZno(k),j)/(Y(CZno(k),CZno(k))+YL(k));
        end
    end
    reduction_mat = [reduction_mat CZno(k)];
    Y = Y1;
end
Y_new = Y;
I_new = I;


% for i = 1:nnew
%     for j = 1:nnew
%         Y_new(i,j) = Y(newno(i),newno(j));
%         % 恒电流节点
%         for k = 1:nCI
%             Y_new(i,j)=Y_new(i,j)-Y(newno(i),CIno(k))*Y(CIno(k),newno(j))/Y(CIno(k),CIno(k));
%         end
%         % 恒阻抗节点
%         for k = 1:nCZ
%             Y_new(i,j)=Y_new(i,j)-Y(newno(i),nCZ(k))*Y(nCZ(k),newno(j))/(Y(nCZ(k),nCZ(k))+YL(k));
%         end
%     end
% end
end

