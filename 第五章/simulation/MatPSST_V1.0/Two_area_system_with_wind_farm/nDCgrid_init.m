function DCgridinit = nDCgrid_init(resultsdc,wB)

busdc=resultsdc.busdc;
branchdc=resultsdc.branchdc;

 %% �γ����Ӿ���T
 
 
    n=size(busdc,1);
    b=size(branchdc,1);
 
 
    TDC=zeros(n,b);%����ڵ�����+1���ѵؽڵ㿼�ǽ�ȥ
    
    
    b_NO   =  linspace(1,b,b).';          %�õ���n����·�ı��
    fdcbus =  branchdc{:,1};%�õ���n����·��Ӧ��frombus�Ľڵ��
    tdcbus =  branchdc{:,2};%�õ���n����·��Ӧ��tobus�Ľڵ��
    
    
    TDC(sub2ind(size(TDC), fdcbus, b_NO))= 1;  %�Ե�b����·��frombus��ֵ1
    TDC(sub2ind(size(TDC), tdcbus, b_NO))=-1;  %�Ե�b����·��tobus��ֵ1
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




