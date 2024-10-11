function [sys,x0,str,ts,simStateCompliance] = TS_net_1_new(t,x,u,flag,busdata,linedata,VSS0,LINE,fbus,ng,Currentin_idx)


persistent YY0 Yff;
persistent YGB YGBd YGBq;

    nBus=size(busdata,1);
    ng2=2*ng;
    ii=1:ng;
    ii2=2*ii;
    ii6=6*ii;
    


%
% The following outlines the general structure of an S-function.
%
switch flag
   
  %%%%%%%%%%%%%%%%%%
  % Initialization %
  %%%%%%%%%%%%%%%%%%
  case 0
    [sys,x0,str,ts,simStateCompliance,YY0,Yff,YGB,YGBd,YGBq]=mdlInitializeSizes(busdata,linedata,VSS0,LINE,nBus,Currentin_idx,fbus,ng);

  %%%%%%%%%%%%%%%
  % Derivatives %
  %%%%%%%%%%%%%%%
  case 1
    sys=[];

  %%%%%%%%%%
  % Update %
  %%%%%%%%%%
  case 2
    sys=[];

  %%%%%%%%%%%
  % Outputs %
  %%%%%%%%%%%
  case 3
    sys=mdlOutputs(t,x,u,YY0,Yff,fbus,ng2,ii,ii2,ii6,YGB,YGBd,YGBq);

  %%%%%%%%%%%%%%%%%%%%%%%
  % GetTimeOfNextVarHit %
  %%%%%%%%%%%%%%%%%%%%%%%
  case 4
    sys=mdlGetTimeOfNextVarHit(t,x,u);

  %%%%%%%%%%%%%
  % Terminate %
  %%%%%%%%%%%%%
  case 9
    sys=[];

  %%%%%%%%%%%%%%%%%%%%
  % Unexpected flags %
  %%%%%%%%%%%%%%%%%%%%
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
function [sys,x0,str,ts,simStateCompliance,YY0,Yff,YGB,YGBd,YGBq]=mdlInitializeSizes(busdata,linedata,VSS0,LINE,nBus,Currentin_idx,fbus,ng)

%
% call simsizes for a sizes structure, fill it in and convert it to a
% sizes array.
%
% Note that in this example, the values are hard coded.  This is not a
% recommended practice as the characteristics of the block are typically
% defined by the S-function parameters.
%
sizes = simsizes;

sizes.NumContStates  = 0;
sizes.NumDiscStates  = 0;
sizes.NumOutputs     = 2*ng;
sizes.NumInputs      = 6*ng+1;
sizes.DirFeedthrough = 1;
sizes.NumSampleTimes = 1;   % at least one sample time is needed

sys = simsizes(sizes);

%
% initialize the initial conditions
%
x0  = [];

%
% str is always an empty matrix
%
str = [];

%
% initialize the array of sample times
%
ts  = [-1 0];

% Specify the block simStateCompliance. The allowed values are:
%    'UnknownSimState', < The default setting; warn and assume DefaultSimState
%    'DefaultSimState', < Same sim state as a built-in block
%    'HasNoSimState',   < No sim state
%    'DisallowSimState' < Error out when saving or restoring the model sim state
simStateCompliance = 'UnknownSimState';

%这里开始计算初始的Y阵
    bus     = busdata;
    branch  = linedata;
    baseMVA = 100;
%% constants
nb = size(bus, 1);          %% number of buses
nl = size(branch, 1);       %% number of lines

%% for each branch, compute the elements of the branch admittance matrix where
%%
%%      | If |   | Yff  Yft |   | Vf |
%%      |    | = |          | * |    |
%%      | It |   | Ytf  Ytt |   | Vt |
%%
stat = branch{:, 'Brc_status'};                    %% ones at in-service branches
Ys = stat ./ (branch{:, 'R'} + 1j * branch{:, 'X'});  %% series admittance
Bc = stat .* branch{:, 'B'};                           %% line charging susceptance
tap = ones(nl, 1);                              %% default tap ratio = 1
i = find(branch{:, 'Tap'});                       %% indices of non-zero tap ratios
tap(i) = branch{i, 'Tap'};                        %% assign non-zero tap ratios
% tap = tap .* exp(1j*pi/180 * branch(:, SHIFT)); %% add phase shifters
Ytt = Ys + 1j*Bc/2;
Yff = Ytt ./ (tap .* conj(tap));
Yft = - Ys ./ conj(tap);
Ytf = - Ys ./ tap;

%% compute shunt admittance
%% if Psh is the real power consumed by the shunt at V = 1.0 p.u.
%% and Qsh is the reactive power injected by the shunt at V = 1.0 p.u.
%% then Psh - j Qsh = V * conj(Ysh * V) = conj(Ysh) = Gs - j Bs,
%% i.e. Ysh = Psh + j Qsh, so ...
Ysh = (bus{:, 'GS'} + 1j * bus{:, 'BS'}) / baseMVA; %% vector of shunt admittances

%% bus indices
f = branch{:, 'F_Bus'};                           %% list of "from" buses
t = branch{:, 'T_Bus'};                           %% list of "to" buses

%% for best performance, choose method based on MATLAB vs Octave and size
if nb < 300 || have_fcn('octave')   %% small case OR running on Octave
    %% build Yf and Yt such that Yf * V is the vector of complex branch currents injected
    %% at each branch's "from" bus, and Yt is the same for the "to" bus end
    i = [1:nl 1:nl]';                           %% double set of row indices
    Yf = sparse(i, [f; t], [Yff; Yft], nl, nb);
    Yt = sparse(i, [f; t], [Ytf; Ytt], nl, nb);

    %% build Ybus
    Ybus = sparse([f;f;t;t], [f;t;f;t], [Yff;Yft;Ytf;Ytt], nb, nb) + ... %% branch admittances
            sparse(1:nb, 1:nb, Ysh, nb, nb);        %% shunt admittance
else                                %% large case running on MATLAB
    %% build connection matrices
    Cf = sparse(1:nl, f, ones(nl, 1), nl, nb);      %% connection matrix for line & from buses
    Ct = sparse(1:nl, t, ones(nl, 1), nl, nb);      %% connection matrix for line & to buses

    %% build Yf and Yt such that Yf * V is the vector of complex branch currents injected
    %% at each branch's "from" bus, and Yt is the same for the "to" bus end
    Yf = sparse(1:nl, 1:nl, Yff, nl, nl) * Cf + sparse(1:nl, 1:nl, Yft, nl, nl) * Ct;
    Yt = sparse(1:nl, 1:nl, Ytf, nl, nl) * Cf + sparse(1:nl, 1:nl, Ytt, nl, nl) * Ct;

    %% build Ybus
    Ybus = Cf' * Yf + Ct' * Yt + ...            %% branch admittances
            sparse(1:nb, 1:nb, Ysh, nb, nb);    %% shunt admittance
end

    nl=linedata{:,'F_Bus'};%frombus
    nr=linedata{:,'T_Bus'};%tobus

    nbr=length(nl);
    nbus=max(max(nl),max(nr));
     
    Ybus=sparse(Ybus);
    
    YY=Ybus;
    
    PPL=busdata{:,'PD'};%PPL=sin(wt)
    QQL=busdata{:,'QD'};
    VVL=VSS0{:,'Vm'};
    YYL=(PPL-1j*QQL)./(VVL.*VVL);

    YY=sparse(YY+diag(YYL));
    
    YY0=YY;
    
    
%%     


    
    R=LINE{1,'R'};
    X=LINE{1,'X'};
    B=LINE{1,'B'}/2;

    F=LINE{1,'F_Bus'};
    T=LINE{1,'T_Bus'};
    

    Yff=zeros(nbus,nbus);
    Yff(F,F)=1/(R+1j*X)+1j*B;
    Yff(F,T)=-1/(R+1j*X);
    Yff(T,F)=-1/(R+1j*X);
    Yff(T,T)=1/(R+1j*X)+1j*B;
    Yff=sparse(Yff);
    

    ng2=2*ng;
    ii=1:ng;
    ii2=2*ii;
%% 
    Y=YY0;
    
    
    
%%         
    Nidx=1:nBus;
    Nidx(Currentin_idx)=[];%GidxNet这个是合并了VSC后的发电机节点号，Nidx看意思应该是所以非发电机节点的节点号
    Ygg=Y(Currentin_idx,Currentin_idx);Ygn=Y(Currentin_idx,Nidx);
    Yng=Y(Nidx,Currentin_idx);Ynn=Y(Nidx,Nidx);
    Yred=Ygg-Ygn/Ynn*Yng;
    
    YGB=zeros(ng2,ng2);
    YGB(ii2-1,ii2-1)=real(Yred(ii,ii));
    YGB(ii2-1,ii2)  =-imag(Yred(ii,ii));
    YGB(ii2,ii2-1)  =imag(Yred(ii,ii));
    YGB(ii2,ii2)    =real(Yred(ii,ii));%
%% 
    YGBd=YGB;
    YGBq=YGB;
if fbus~=999
    Y(fbus,fbus)=Y(fbus,fbus)+1;
%      Y(8,8)=Y(8,8)+1e5;
    
    Nidx=1:nBus;
    Nidx(Currentin_idx)=[];
    Ygg=Y(Currentin_idx,Currentin_idx);Ygn=Y(Currentin_idx,Nidx);
    Yng=Y(Nidx,Currentin_idx);Ynn=Y(Nidx,Nidx);
    Yred=Ygg-Ygn/Ynn*Yng;
    
    YGBd=zeros(ng2,ng2);
    YGBd(ii2-1,ii2-1)=real(Yred(ii,ii));
    YGBd(ii2-1,ii2)  =-imag(Yred(ii,ii));
    YGBd(ii2,ii2-1)  =imag(Yred(ii,ii));
    YGBd(ii2,ii2)    =real(Yred(ii,ii));%
%% 计算切除短路后的矩阵
    Y=Y-Yff;
    Nidx=1:nBus;
    Nidx(Currentin_idx)=[];
    Ygg=Y(Currentin_idx,Currentin_idx);Ygn=Y(Currentin_idx,Nidx);
    Yng=Y(Nidx,Currentin_idx);Ynn=Y(Nidx,Nidx);
    Yred=Ygg-Ygn/Ynn*Yng;
    
    YGBq=zeros(ng2,ng2);
    YGBq(ii2-1,ii2-1)=real(Yred(ii,ii));
    YGBq(ii2-1,ii2)  =-imag(Yred(ii,ii));
    YGBq(ii2,ii2-1)  =imag(Yred(ii,ii));
    YGBq(ii2,ii2)    =real(Yred(ii,ii));%
end

% end mdlInitializeSizes

%
%=============================================================================
% mdlOutputs
% Return the block outputs.
%=============================================================================
%
function sys=mdlOutputs(t,x,u,YY0,Yff,fbus,ng2,ii,ii2,ii6,YGB,YGBd,YGBq)

Y=YY0;
%%            =========Fault time===========

t1=5;    %5
t2=50000;    

t3=50000;  

%%%% 

if (u(1)>=t1)&&(u(1)<t2)
    
    YGB=YGBd;
    
end

if (u(1)>=t2)&&(u(1)<t3) 
    
    YGB=YGBq;
         
end



%%%%

YGB(  (ng2+1).*(ii2-1)-ng2  )      =   YGB(  (ng2+1).*(ii2-1)-ng2  ) +u(ii6-2).';
YGB( (ng2+1).*(ii2)-ng2-1   )      =   YGB( (ng2+1).*(ii2)-ng2-1   ) +u(ii6).';
YGB( (ng2+1).*(ii2-1)-ng2+1 )      =   YGB( (ng2+1).*(ii2-1)-ng2+1 ) +u(ii6+1).';
YGB( (ng2+1).*(ii2)-ng2     )      =   YGB( (ng2+1).*(ii2)-ng2     ) +u(ii6-1).';


%%%%
 Ixy=zeros(ng2,1);
 Ixy([ii2-1,ii2])=u([ii6-4,ii6-3]);


Vxy=YGB\Ixy;


%%%%
% Vg=zeros(ng,1);
% for i=1:ng
%     Vg(i)=Vxy(2*i-1)+1j*Vxy(2*i);
% end

%Vn=-Ynn\Yng*Vg;%
%Vnang=angle(Vn);


sys = [Vxy];

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

