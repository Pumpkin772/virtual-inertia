function [resultsac,resultsdc,converged] = powerflowcaculate...
    (mpcac,mpcdc,PaggWind,AggWind_in_gidx,SynMach_in_gidx)

%% 
    pfopt1 = macdcoption;
    pfopt1(13)=1;
    pfopt1([2,4,6,8])=40;
    pfopt0 = macdcoption;
    pfopt0(13)=0;
    pfopt0([2,4,6,8])=40;
%% 
tablempcac=mpcTVSA(mpcac,'Mac_A2T');
tablempcdc=mpcTVSA(mpcdc,'Mdc_A2T');


%% 
tablempcac.gen{cellstr(string(AggWind_in_gidx)),'PG'} =  PaggWind;

    
%% 
    [resultsac, resultsdc, converged]    =   runacdcpf(mpcTVSA(tablempcac,'AC_T2A'),mpcTVSA(tablempcdc,'DC_T2A'));

    
    resultsac=mpcTVSA(resultsac,'Rac_A2T');
    resultsdc=mpcTVSA(resultsdc,'Rdc_A2T');

%% 

Gen_idx=zeros(size(resultsac.gen,1),1);
GWType=strings(size(resultsac.gen,1),1);

% C = cell(size(resultsac.gen,1),2);
% Ctable=cell2table(C,...
%     'VariableNames',{'Gen_idx','GWType'});

% resultsac.gen=addvars(resultsac.gen,Gen_idx,GWType,'NewVariableNames',);

resultsac.gen{cellstr(string(SynMach_in_gidx)),'Gen_idx'}=[1:size(SynMach_in_gidx,2)].';
resultsac.gen{cellstr(string(SynMach_in_gidx)),'GWType'}=string('Gen');

resultsac.gen{cellstr(string(AggWind_in_gidx)),'Gen_idx'}=[1:size(AggWind_in_gidx,1)].';
resultsac.gen{cellstr(string(AggWind_in_gidx)),'GWType'}=string('AggWind');


Conv_idx=[1:size(resultsdc.convdc,1)].';
% resultsdc.convdc=addvars(resultsdc.convdc,Conv_idx,'NewVariableNames','Conv_idx');
resultsdc.convdc{cellstr(string(Conv_idx)),'Conv_idx'}=Conv_idx;
end



















