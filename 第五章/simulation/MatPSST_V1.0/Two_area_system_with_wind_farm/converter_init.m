function [VSCinit] = converter_init(resultsac,resultsdc)
CP       =  1;%control_P 
CUDC     =  2;%control_Udc
CQ       =  1;%control_Q
CUAC     =  2;%control_Uac
        
w        =  1;
% 
% d: 1==CP 2==CUDC
% q: 1==CQ 2==CUAC
d_axis_control=[CP,CUDC,CP];
q_axis_control=[CQ,CQ,CQ];

if any(((d_axis_control==3)&(q_axis_control==3))~=((d_axis_control==3)|(q_axis_control==3)))
    error('auto模式dq必须一起选')
end


VSCinit.d_axis_control=d_axis_control;
VSCinit.q_axis_control=q_axis_control;
%% 
%  1VSC序号   2 Vdcbase     3 Vacbase        4            5     6 SvscBi
VSC_data=...
[
     1           400          200           0.001     0.0002       100;
     2           400          200           0.001     0.0002       100;
     3           400          200           0.001     0.0002       100; 

];

Vdcbase    =          VSC_data(:,2);
Vacbase    =          VSC_data(:,3);

ConvandBusdc=join(resultsdc.convdc,resultsdc.busdc,'LeftKeys','Conv_Bus','RightKeys','Busdc_I');
% con_dcidx_con          =     resultsdc.convdc(:,1);
% [~,~,con_dchidx_bus]   =     intersect(con_dcidx_con,resultsdc.busdc(:,1));

Lpr        =          ConvandBusdc{:,'Xconv'};
Rpr        =          ConvandBusdc{:,'Rconv'};
CeqA       =          ConvandBusdc{:,'Cdc'};


SvscBi     =          VSC_data(:,6);

VSCinit.Lpr    =    Lpr;
VSCinit.Rpr    =    Rpr;
VSCinit.CeqA   =    CeqA;
VSCinit.Vdcbase=    Vdcbase;
VSCinit.Vacbase=    Vacbase;
VSCinit.SvscBi =    SvscBi;
%% 
%               23
%    VSC序号   PLLPI
VSC_PLL_data=...
    [
    1          100,6000
    2          100,6000
    3          100,6000
    ];

VSCinit.Kppll=VSC_PLL_data(:,2);
VSCinit.Kipll=VSC_PLL_data(:,3);

%% 
%       1         2             34         56          78         9,10       11,12        13,14        
%    VSC序号     PWMTo       d轴内环PI   q轴内环PI    d外VdcPI    d外PPI      q外VacPI      q外QPI      
VSC_Control_data=...
    [
    1       0.001       0.04,0.28   0.04,0.28    45,950     5.06,382.07     10,200    -5.06,-150.00      
    2       0.001       0.04,0.28   0.04,0.28    45,950     5.06,382.07     10,200    -5.06,-150.00      
    3       0.001       0.04,0.28   0.04,0.28    45,950     5.06,382.07     10,200    -5.06,-150.00      
    ];
% 4       0.001       0.04,0.28   0.04,0.28    25,950     5.06,382.07     10,200    -5.06,-382.07    20,100     10,10;
%
VSCinit.To     =          VSC_Control_data(:,2);

VSCinit.KpId   =          VSC_Control_data(:,3);
VSCinit.KiId   =          VSC_Control_data(:,4);

VSCinit.KpIq   =          VSC_Control_data(:,5);
VSCinit.KiIq   =          VSC_Control_data(:,6);



for ii=1:1:size(VSC_Control_data,1)
    switch d_axis_control(ii)
        case CP
        case CUDC
            VSCinit.KpId(ii)      =          0.1;
            VSCinit.KiId(ii)      =          1;
            
            VSCinit.KpIq(ii)      =          0.1;
            VSCinit.KiIq(ii)      =          1;
        case AUTO
            VSCinit.KpId(ii)      =          2;
            VSCinit.KiId(ii)      =          1;
            
            VSCinit.KpIq(ii)      =          2;
            VSCinit.KiIq(ii)      =          1;
        otherwise
            error('invaild control strategy');
    end
end
%  d-axis control 
%
VSCinit.Kpdc =        VSC_Control_data(:,7);
VSCinit.Kidc =        VSC_Control_data(:,8);

%
VSCinit.KpP  =        VSC_Control_data(:,9);
VSCinit.KiP  =        VSC_Control_data(:,10);


%  q-axis control 
%
VSCinit.Kpac =        VSC_Control_data(:,11);
VSCinit.Kiac =        VSC_Control_data(:,12);
% 
VSCinit.KpQ  =        VSC_Control_data(:,13);
VSCinit.KiQ  =        VSC_Control_data(:,14);






%%

VSCinit.Qref        =        -resultsdc.convdc{:,5}./SvscBi;
VSCinit.Pref        =        -resultsdc.convdc{:,4}./SvscBi;
VSCinit.Uref        =         resultsdc.convdc{:,30};

VSCinit.Udcref      =         ConvandBusdc{:,'Vdc'};

%%
        UeqA       =        ConvandBusdc{:,'Vdc'};

VSCinit.UeqA       =        UeqA;%
%% 

PCC_idx            =    ConvandBusdc{:,'Busac_I'};%
Angle_PCC          =    resultsac.bus{cellstr(string(PCC_idx)),9}*pi/180;




S_inject_insystem      =    -(resultsdc.convdc{:,4}+1j*resultsdc.convdc{:,5})./SvscBi;
Amplitude_PCC          =    resultsac.bus{cellstr(string(PCC_idx)),8};

Usd0                =    Amplitude_PCC;   
Usq0                =    zeros(size(resultsdc.convdc,1),1);
Usdq0               =    Usd0+1j.*Usq0;   

Isdq0              =    conj(S_inject_insystem./Usdq0);
Isd0                =    real(Isdq0);
Isq0                =    imag(Isdq0);


S_inject_insystem      =    -(resultsdc.convdc{:,4}+1j*resultsdc.convdc{:,5})./SvscBi;
Voltage_Phasor_PCC     =    resultsac.bus{cellstr(string(PCC_idx)),8}.*exp(1j*Angle_PCC);
Current_Phasor_to_PCC  =    conj(S_inject_insystem./Voltage_Phasor_PCC);

	    Ucd0        =    Usd0-Rpr.*Isd0+Lpr.*Isq0;
        Ucq0        =    Usq0-Rpr.*Isq0-Lpr.*Isd0;
        
        Uc          =    sqrt(Ucd0.^2+Ucq0.^2);
        
    



VSCinit.Dels_hvdc      =    Angle_PCC;
VSCinit.Isx0        =    real(Current_Phasor_to_PCC);
VSCinit.Isy0        =    imag(Current_Phasor_to_PCC);

VSCinit.Ucd0        =    Ucd0;
VSCinit.Ucq0        =    Ucq0;

VSCinit.Md0         =    Usd0-Ucd0+w.*Lpr.*Isq0;
VSCinit.Mq0         =    Usq0-Ucq0-w.*Lpr.*Isd0;
VSCinit.Nd0         =    Isd0;
VSCinit.Nq0         =    Isq0;
%% end

end








