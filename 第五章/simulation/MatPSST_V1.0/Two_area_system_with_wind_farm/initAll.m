%% Common parameters

w         =          1;
fB        =          50;                   
wB        =          2*pi*fB;



SynMach_in_gidx     =    [1 2 3 4 12];%The node number of the synchronous motor connected to the large power grid. 
                                    % Please arrange it in the order in the gen matrix.   
AggWind_in_gidx     =    [8];%The node number of the wind turbine connected to the large power grid. Please arrange it in the order in the gen matrix.          

%  Single fan capacity 
AggSwgB           =    [2];%2MVA
%  Number of aggregate fans 
AggNN             =    [50];%a DFIG model represent 50 DFIGs
%  Rotor Speed of DFIG 
AggWr0            =    [0.9];%
%  Initial power of fan 
AggPwind_S   =    AggWr0.*...
    1.162.*AggWr0.^2;
AggPwind_A   =    AggPwind_S.*AggNN.*AggSwgB;






mpcdc =loadcasedc('P_DCdata');
mpcac =loadcase('P_ACdata');




[resultsac,resultsdc,converged]=powerflowcaculate(mpcac,mpcdc,AggPwind_A,AggWind_in_gidx,SynMach_in_gidx);

  

if(converged~=1)
    error('潮流计算不收敛')
end
clear converged

SgridB    =     resultsac.baseMVA;
SdcgridB  =     resultsdc.baseMVAdc;




%% AC network initialization 

fbranchac = resultsac.branch(1,:);  % Set up short line 

[ACgirdinit,opdata]=ACgrid_init(resultsac,resultsdc,fbranchac);
 

%% DC network parameter initialization 
DCgridinit=nDCgrid_init(resultsdc,wB);
%% Converter parameter initialization 
VSCinit=converter_init(resultsac,resultsdc);
% 
%% DFIG initialization 
[Windinit,AVm0] = nWind_machine_init(opdata,AggWind_in_gidx,AggWr0,AggNN,AggSwgB);


%% Synchronous motor initialization 
Geninit=Generator_init(opdata,SynMach_in_gidx);


