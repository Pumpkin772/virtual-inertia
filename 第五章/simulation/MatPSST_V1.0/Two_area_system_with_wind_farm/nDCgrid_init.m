function DCgridinit = nDCgrid_init(resultsdc,wB)

busdc=resultsdc.busdc;
branchdc=resultsdc.branchdc;

 %% 形成连接矩阵T
 
 
    n=size(busdc,1);
    b=size(branchdc,1);
 
 
    TDC=zeros(n,b);%这里节点数量+1，把地节点考虑进去
    
    
    b_NO   =  linspace(1,b,b).';          %得到第n条线路的编号
    fdcbus =  branchdc{:,1};%得到第n条线路对应的frombus的节点号
    tdcbus =  branchdc{:,2};%得到第n条线路对应的tobus的节点号
    
    
    TDC(sub2ind(size(TDC), fdcbus, b_NO))= 1;  %对第b条线路的frombus赋值1
    TDC(sub2ind(size(TDC), tdcbus, b_NO))=-1;  %对第b条线路的tobus赋值1
    TDCt   =  TDC.';


DCgridinit.TDCt=TDCt;
DCgridinit.TDC=TDC;

%% 
 line_resistor_vector       =        resultsdc.branchdc{:,3};
 line_inductor_vector       =        resultsdc.branchdc{:,4};
 


               
DCgridinit.RdcA=line_resistor_vector;
DCgridinit.LdcA=line_resistor_vector;
%% 

    linepower_frombus       =        resultsdc.branchdc{:,10}/resultsdc.baseMVAdc;
    frombus_idx             =        resultsdc.branchdc{:,1};
    busDCvoltagevalue       =        resultsdc.busdc{frombus_idx,5};
          linecurrent       =        linepower_frombus./busDCvoltagevalue;
          



DCgridinit.IdcA=linecurrent;
%% 

busDC_Capacitor_value       =        resultsdc.busdc(:,9);

                 
DCgridinit.CdcA=busDC_Capacitor_value;                        
%%                  
    busDCvoltagevalue       =        resultsdc.busdc(:,5);
    
    
 DCgridinit.UdcA=busDCvoltagevalue;
 
end




