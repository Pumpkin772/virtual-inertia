function [resultsac,resultsdc,converged] = powerflowcaculate...
    (mpcac,mpcdc)

%% 
    pfopt1 = macdcoption;
    pfopt1(31)=0;
%     pfopt1([2,4,6,8])=40;
%     pfopt0 = macdcoption;
%     pfopt0(13)=0;
%     pfopt0([2,4,6,8])=40;
%% 
tablempcac=mpcTVSA(mpcac,'Mac_A2T');
tablempcdc=mpcTVSA(mpcdc,'Mdc_A2T');


    
%% 
    [resultsac, resultsdc, converged]    =   runacdcpf(mpcTVSA(tablempcac,'AC_T2A'),mpcTVSA(tablempcdc,'DC_T2A'));

    
    resultsac=mpcTVSA(resultsac,'Rac_A2T');
    resultsdc=mpcTVSA(resultsdc,'Rdc_A2T');

%% 



end



















