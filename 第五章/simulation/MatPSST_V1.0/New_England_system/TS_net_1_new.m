function [sys,x0,str,ts,simStateCompliance] = TS_net_1_new(t,x,u,flag,YY_load,ng)


persistent Yff 

m=29;
n=39;



%
% The following outlines the general structure of an S-function.
%
switch flag
  case 0
    [sys,x0,str,ts,simStateCompliance,Yff]=mdlInitializeSizes(YY_load,m,n,ng);
  case 1
    sys=[];
  case 2
    sys=[];
  case 3
    sys=mdlOutputs(t,x,u,Yff,YY_load,m,n);
  case 4
    sys=mdlGetTimeOfNextVarHit(t,x,u);  
  case 9
    sys=[];
  otherwise
    DAStudio.error('Simulink:blocks:unhandledFlag', num2str(flag));
end

% end sfuntmpl











%
%=============================================================================
% mdlInitializeSizes
% Return the sizes, initial conditions, and sample times for the S-function.
%=============================================================================
%
function [sys,x0,str,ts,simStateCompliance,Yff]=mdlInitializeSizes(YY_load,m,n,ng)


sizes = simsizes;

sizes.NumContStates  = 0;
sizes.NumDiscStates  = 0;
sizes.NumOutputs     = 2*ng+39-ng;
sizes.NumInputs      = 6*ng+1;
sizes.DirFeedthrough = 1;
sizes.NumSampleTimes = 1;

sys = simsizes(sizes);


x0  = [];


str = [];

ts  = [-1 0];

simStateCompliance = 'UnknownSimState';



    R=0.00130;
    X=0.02130;
    B=0.22140/2;
    
    Yff=zeros(n,n);
    Yff(3,3)=1/(R+1i*X)+1i*B;
    Yff(3,4)=-1/(R+1i*X);
    Yff(4,3)=-1/(R+1i*X);
    Yff(4,4)=1/(R+1i*X)+1i*B;
    
    


function sys=mdlOutputs(t,x,u,Yff,YY_load,m,n)
Y=YY_load;
t=u(1);
% ===============Fault time==================
t1=5;    
t2=5.1;    
t3=5.2;    

% ==============================================
%%%% 
if (t>=t1)&&(t<t2)
  Y(3,3)=Y(3,3)+1e5;  
end

%%%% 
if (t>=t2)&&(t<t3)     
    Y=Y-Yff;   
end

%%%%%%  
Ygg=Y(m+1:n,m+1:n);
Ygn=Y(m+1:n,1:m);
Yng=Y(1:m,m+1:n);
Ynn=Y(1:m,1:m);

Yred=Ygg-Ygn*inv(Ynn)*Yng;

GidxNet=[30:39];
%% 
Nidx=1:n;Nidx(GidxNet)=[];
Ygg=Y(GidxNet,GidxNet);Ygn=Y(GidxNet,Nidx);
Yng=Y(Nidx,GidxNet);Ynn=Y(Nidx,Nidx);
Yred=Ygg-Ygn/Ynn*Yng;


%% 
ng=length(Yred);

YGB=zeros(2*ng,2*ng);

for i=1:ng
    for k=1:ng
        YGB(2*i-1,2*k-1)=real(Yred(i,k));
        YGB(2*i-1,2*k)=-imag(Yred(i,k));
        YGB(2*i,2*k-1)=imag(Yred(i,k));
        YGB(2*i,2*k)=real(Yred(i,k));
    end
end

%% 
for i=1:ng
    YGB(2*i-1,2*i-1)=YGB(2*i-1,2*i-1)+u(6*i-2);
    YGB(2*i-1,2*i)=YGB(2*i-1,2*i)+u(6*i);
    YGB(2*i,2*i-1)=YGB(2*i,2*i-1)+u(6*i+1);
    YGB(2*i,2*i)=YGB(2*i,2*i)+u(6*i-1);
end

%% 
Ixy=[];
for i=1:ng
    Ixy=[Ixy;u(6*i-4);u(6*i-3)];
end

Vxy=YGB\Ixy;

%% 
Vg=zeros(ng,1);
for i=1:ng
    Vg(i)=Vxy(2*i-1)+j*Vxy(2*i);
end

Vn=Ynn\Yng*Vg;

%% 
sys=[Vxy;abs(Vn)];


% end mdlOutputs

%
%=============================================================================
% mdlGetTimeOfNextVarHit
% Return the time of the next hit for this block.  Note that the result is
% absolute time.  Note that this function is only used when you specify a
% variable discrete-time sample time [-2 0] in the sample time array in
% mdlInitializeSizes.
%=============================================================================
%
function sys=mdlGetTimeOfNextVarHit(t,x,u)

sampleTime = 1;    %  Example, set the next hit to be one second later.
sys = t + sampleTime;

% end mdlGetTimeOfNextVarHit

