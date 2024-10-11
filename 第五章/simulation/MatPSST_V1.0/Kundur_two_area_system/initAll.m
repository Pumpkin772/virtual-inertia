%% Common parameters

w         =          1;
fB        =          50;                   
wB        =          2*pi*fB;



SynMach_in_gidx     =    [1 2 3 4];%The node number of the synchronous motor connected to the large power grid. 
                                    % Please arrange it in the order in the gen matrix.   







mpcdc =loadcasedc('P_DCdata');
mpcac =loadcase('P_ACdata');




[resultsac,resultsdc,converged]=powerflowcaculate(mpcac,mpcdc,SynMach_in_gidx);

  

if(converged~=1)
    error('潮流计算不收敛')
end
clear converged

SgridB    =     resultsac.baseMVA;
SdcgridB  =     resultsdc.baseMVAdc;




%% AC network initialization 

fbranchac = resultsac.branch(1,:);  % Set up short line 

[ACgirdinit,opdata]=ACgrid_init(resultsac,resultsdc,fbranchac);
 



%% Synchronous motor initialization 
Geninit=Generator_init(opdata,SynMach_in_gidx);


