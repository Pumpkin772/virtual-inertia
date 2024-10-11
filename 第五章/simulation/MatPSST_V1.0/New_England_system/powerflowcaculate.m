%*********************************************************************
%*** New England 10 machine 39 bus 46 line power system *************
%*********************************************************************
function [opdata,YY_load]=powerflowcaculate(busdata,linedata)

% global YY0 PF opdata YY_load; %Y0Ϊ��������������ɵ�Y��

%  ����
j=sqrt(-1);


%�����������
lta=1e-5;

%%%  ��⵼�ɾ���
nl=linedata(:,1);
nr=linedata(:,2);
R=linedata(:,3);
X=linedata(:,4);
B=linedata(:,5)/2; %%ע�⣬linedata�и�����BΪ��·�ֲܷ�����
Tap=linedata(:,6);

nbr=length(nl);
nbus=max(max(nl),max(nr));
Z=R+j*X;
y=ones(nbr,1)./Z; %%֧·�迹ȡ����

Y=zeros(nbus,nbus);

for k=1:nbr
    if (nl(k)>0)&(nr(k)>0)
        if Tap(k)==0
            Y(nl(k),nr(k))=Y(nl(k),nr(k))-y(k);
        else
            Y(nl(k),nr(k))=Y(nl(k),nr(k))-y(k)/Tap(k);
        end
        Y(nr(k),nl(k))=Y(nl(k),nr(k));
    end
end

for i=1:nbus
    for k=1:nbr
        if nl(k)==i|nr(k)==i
            if (Tap(k)~=0)&(nl(k)==i)
                Y(i,i)=Y(i,i)+y(k)/(Tap(k)*Tap(k))+j*B(k);
            else
                Y(i,i)=Y(i,i)+y(k)+j*B(k);
            end
        end
    end
end

Gs=busdata(:,8);
Bs=busdata(:,9);
Ys=Gs+j*Bs;
for i=1:nbus
    Y(busdata(i,1),busdata(i,1))=Y(busdata(i,1),busdata(i,1))+Ys(i);
end

Y=sparse(Y);

%%%%******************************************
%%% m��PQ�ڵ㣬n-m-1��PV�ڵ㣬n��Ϊƽ��ڵ�
%%%%******************************************
Busno=busdata(:,10);
PQno=find(busdata(:,10)==3);    
PVno=find(busdata(:,10)==2);
SWno=find(busdata(:,10)==1);
PQVno=[PQno;PVno];
PQVSno=[PQno;PVno;SWno];

m=length(PQno);
n=length(Busno);

%%%���սڵ�����ڵ㵼�ɾ���
%%%***sy��Ϊ10��39�ڵ�����ڵ�˳����ǰ���PQ,PV,SW���ģ������ⲽûӰ��
YY0=Y;
Y=[YY0(busdata(PQno,1),:);YY0(busdata(PVno,1),:);YY0(busdata(SWno,1),:)];
YY0=Y;
Y=[YY0(:,busdata(PQno,1)) YY0(:,busdata(PVno,1)) YY0(:,busdata(SWno,1))];
YG=real(Y);
YB=imag(Y);


% %��ʼ���� %%***sy�ڵ㹦����ע��ڵ�Ϊ��
Pg=busdata(PQVno,4);
PL=busdata(PQVno,6);
Ps=Pg-PL;  

Qg=busdata(PQno,5);
QL=busdata(PQno,7);
Qs=Qg-QL;

Vs=busdata(PQVSno,2);
Dels=busdata(PQVSno,3)*pi/180;

P=zeros(n-1,1);
Q=zeros(m,1);

for i=1:n-1
    P(i)=0;
    for k=1:n
        P(i)=P(i)+Vs(i)*Vs(k)*(YG(i,k)*cos(Dels(i)-Dels(k))+YB(i,k)*sin(Dels(i)-Dels(k)));
    end
end

for i=1:m
    Q(i)=0;
    for k=1:n
        Q(i)=Q(i)+Vs(i)*Vs(k)*(YG(i,k)*sin(Dels(i)-Dels(k))-YB(i,k)*cos(Dels(i)-Dels(k)));
    end
end

dPs=Ps-P;
dQs=Qs-Q;
dPQ=[dPs;dQs];

dPQmax=max(abs(dPQ));

flag=0;

%%%ţ�ٷ�������⳱��
while (dPQmax>lta)
    
    Jac=zeros(n-1+m);
    JH=zeros(n-1);
    JN=zeros(n-1,m);
    JK=zeros(m,n-1);
    JL=zeros(m,m);
    
    %%% ����JH
    for i=1:n-1
        for k=1:n-1
            if i==k
                QQs=0;
                for l=1:n
                    QQs=QQs+Vs(i)*Vs(l)*(YG(i,l)*sin(Dels(i)-Dels(l))-YB(i,l)*cos(Dels(i)-Dels(l)));
                end
                JH(i,k)=(Vs(i)^2)*YB(i,k)+QQs;
            else
                JH(i,k)=-Vs(i)*Vs(k)*(YG(i,k)*sin(Dels(i)-Dels(k))-YB(i,k)*cos(Dels(i)-Dels(k)));
            end
        end
    end
    
    %%% ����JN
    for i=1:n-1
        for k=1:m
            if i==k
                PPs=0;
                for l=1:n
                    PPs=PPs+Vs(i)*Vs(l)*(YG(i,l)*cos(Dels(i)-Dels(l))+YB(i,l)*sin(Dels(i)-Dels(l)));
                end
                JN(i,k)=-(Vs(i)^2)*YG(i,k)-PPs;
            else
                JN(i,k)=-Vs(i)*Vs(k)*(YG(i,k)*cos(Dels(i)-Dels(k))+YB(i,k)*sin(Dels(i)-Dels(k)));
            end
        end
    end
    
    %%% ����JK
    for i=1:m
        for k=1:n-1
            if i==k
                PPs=0;
                for l=1:n
                    PPs=PPs+Vs(i)*Vs(l)*(YG(i,l)*cos(Dels(i)-Dels(l))+YB(i,l)*sin(Dels(i)-Dels(l)));
                end
                JK(i,k)=(Vs(i)^2)*YG(i,k)-PPs;
            else
                JK(i,k)=Vs(i)*Vs(k)*(YG(i,k)*cos(Dels(i)-Dels(k))+YB(i,k)*sin(Dels(i)-Dels(k)));
            end
        end
    end
    
    %%% ����JL
    for i=1:m
        for k=1:m
            if i==k
                QQs=0;
                for l=1:n
                    QQs=QQs+Vs(i)*Vs(l)*(YG(i,l)*sin(Dels(i)-Dels(l))-YB(i,l)*cos(Dels(i)-Dels(l)));
                end
                JL(i,k)=(Vs(i)^2)*YB(i,k)-QQs;
            else
                JL(i,k)=-Vs(i)*Vs(k)*(YG(i,k)*sin(Dels(i)-Dels(k))-YB(i,k)*cos(Dels(i)-Dels(k)));
            end
        end
    end
    
    Jac=[JH JN;JK JL];
    dPQ=[dPs;dQs];
    
    dDV=-inv(Jac)*dPQ;
    dD=dDV(1:n-1);
    dVs=diag(Vs(1:m))*dDV(n:n+m-1);
    
    Dels(1:n-1)=Dels(1:n-1)+dD;
    Vs(1:m)=Vs(1:m)+dVs;
    
    %%%%%����dPQs
    for i=1:n-1
        P(i)=0;
        for k=1:n
            P(i)=P(i)+Vs(i)*Vs(k)*(YG(i,k)*cos(Dels(i)-Dels(k))+YB(i,k)*sin(Dels(i)-Dels(k)));
        end
    end
    
    for i=1:m
        Q(i)=0;
        for k=1:n
            Q(i)=Q(i)+Vs(i)*Vs(k)*(YG(i,k)*sin(Dels(i)-Dels(k))-YB(i,k)*cos(Dels(i)-Dels(k)));
        end
    end
    
    dPs=Ps-P;
    dQs=Qs-Q;
    dPQ=[dPs;dQs];
    
    dPQmax=max(abs(dPQ));
    
    VDelta=[busdata(PQVSno,1) Vs Dels];
    
    flag=flag+1;
    if flag>20
        disp(['***********************************'])
        disp(['There is no solution for power flow'])
        disp(['***********************************'])
        break
    end     
end

%����ƽ��ڵ��ע�빦��
i=n;
Qsn=0;
for l=1:n
    Qsn=Qsn+Vs(i)*Vs(l)*(YG(i,l)*sin(Dels(i)-Dels(l))-YB(i,l)*cos(Dels(i)-Dels(l)));
end

Psn=0;
for l=1:n
    Psn=Psn+Vs(i)*Vs(l)*(YG(i,l)*cos(Dels(i)-Dels(l))+YB(i,l)*sin(Dels(i)-Dels(l)));
end

%===============ƽ��ڵ�����ڹ���===============
Sn=Psn+j*Qsn;
Vss=Vs.*cos(Dels)+j*Vs.*sin(Dels);
%=======��PQ,PV,Swing�Ų��Ľڵ��ѹ��Ϣ===========
%VSS=[��ţ���ѹ��ֵ�����(��)��������������]
Vss=Vs.*cos(Dels)+j*Vs.*sin(Dels);
VSS=[busdata(PQVSno,1) abs(Vss)  Dels*180/pi Vss];
% VSS=[busdata(PQVSno,1) abs(Vss)  Dels Vss];

VSS0=[];%% VSS0ÿ������Ϊ�ڵ��ţ���ѹ��ֵ����ǣ�������ʽ�ĵ�ѹ_sy
% for i=1:n
%    for k=1:n
%        if VSS(k,1)==i
%            VSS0=[VSS0;VSS(k,:)];
%        end
%    end
% end
%==============��˳���ŵĽڵ��ѹ��Ϣ=============
%VSS0=[��ţ���ѹ��ֵ�����(��)��������������]
VSS0=sortrows(VSS,1);

%���㷢����ڵ���޹�����
QQg=zeros(length(PVno),1);
for i=1:length(PVno);
    QQg(i)=0;
    for l=1:n
        QQg(i)=QQg(i)+Vs(m+i)*Vs(l)*(YG(m+i,l)*sin(Dels(m+i)-Dels(l))-YB(m+i,l)*cos(Dels(m+i)-Dels(l)));
    end
end

QQG=[busdata(PVno,1) QQg]; %QQg�㵽�����1~9���޹�Q=Qg-QL��sy��

% �����ƽ��ڵ������нڵ��޹�Q=Qg-QL
for ii=1:n-1
    Qs(ii)=0;
    for l=1:n
        Qs(ii)=Qs(ii)+Vs(ii)*Vs(l)*(YG(ii,l)*sin(Dels(ii)-Dels(l))-YB(ii,l)*cos(Dels(ii)-Dels(l)));
    end
end

%=================��˳���ŵĳ�����Ϣ==================
%PF=[��ţ���ѹ��ֵ�����(��)����������������Pg,Qg,PL,QL]
PL=busdata(PQVSno,6);
QL=busdata(PQVSno,7);
%%%(sy)PF�в�����ƽ��ڵ��P,Q������PV�ڵ��Q������Ҫ��Qg��PgҪ���ϸ��ɵġ���ΪPs=Pg-PL;Qs=Qg-QL
Pg=[busdata(1:38,4);Psn+busdata(39,6)];
Qg=[Qs+busdata(1:38,7);Qsn+busdata(39,7)];
PQ=[Pg Qg PL QL];
disp('      ���               ��ѹ��ֵ                ���(��)              ������������           Pg             Qg             PL           QL')
PF=[VSS0 PQ]

PPL=busdata(:,6);
QQL=busdata(:,7);
VVL=VSS0(:,2);
YYL=(PPL-j*QQL)./(VVL.*VVL);
%%% �����ɹ鲢�����ɾ�����ȥ
YY=sparse(YY0+diag(YYL));
YY_load=full(YY);

%=================��˳���ŵķ������Ϣ==================
SG10=busdata(n,6)+j*busdata(n,7)+Sn;
SG=[busdata(PVno,1) busdata(PVno,4) QQg+busdata(PVno,7);n real(SG10) imag(SG10)];
%%%%%%(sy)opdata�в�����ƽ��ڵ��P,Q������PV�ڵ��Q������Ҫ��Qg��PgҪ���ϸ��ɵġ���ΪPs=Pg-PL;Qs=Qg-QL
%opdata=[��ţ�Pg��Qg����ѹ��ֵ�����(��)]
opdata=[SG VSS0(30:39,2:3)]

%=================��˳���ŵ�֧·��Ϣ==================
% Sline1=VSS0(nl,4).*conj((VSS0(nl,4)-VSS0(nr,4)).*y+VSS0(nl,4).*(j*B));
% Sline2=VSS0(nr,4).*conj((VSS0(nr,4)-VSS0(nl,4)).*y+VSS0(nr,4).*(j*B));
% line_Power=[nl nr real(Sline1) imag(Sline1) nr nl real(Sline2) imag(Sline2)];
% Ploss=sum(real(Sline1))+sum(real(Sline2));


end
