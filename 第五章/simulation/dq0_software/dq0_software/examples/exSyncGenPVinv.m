% Dynamics and stability analysis of large power systems
% with physical synchronous machines & a photovoltaic inverter.

% To begin:
% -- choose a network below ('choose_network'),
% -- select which generator will be replaced with a photovoltaic
% inverter (variable 'synrp')

% Note: the photovoltaic inverter rated power and input power may be adjusted
% using the variables 'Pinv_rt' and 'Pdc_op'.

% Short description:
% a) The system is connected to an infinite bus.
% b) physical model of synchronous machines, as described in "Electric
% Machinery" by Fitzgerald.
% c) all DQ0 signals are referenced to infinite bus
% d) photovoltaic inverter: typical dq0-based PQ control scheme similar to the one
% described in N. Kroutikova, C. A. Hernandez-Aramburo, and T. C. Green,
% “Statespace model of grid-connected inverters under current control mode,”
% IET Electr. Power App., vol. 1, no. 3, pp. 329–338, 2007.
%
% ** Associated Simulink file: 'exSyncGenPVinv_sim'
%
% the MatPower package should be installed (available online).

clc;
close all;

choose_network = 1;
% Choose a network from the list below:

synrp = 2; 
% Choose which synchronous machine to replace with a
% photovoltaic inverter. (This is the machine number, not the bus number)
% for example, in the 30-bus system, bus 1 is the
% infinite bus, and buses 2,5,8,11,13 represent synchronous generators.
% In this example 'synrp' must be in the range 1-5. If for instance
% synrp=4, then the synchronous machine on bus 11 is replaced by a photovoltaic
% inveter.

define_constants; % MatPower function. Defines MatPower constants.

%%%%%%%%%% Simulation parameters %%%%%%%%%%%%%
sim_time = 20;  % [sec] simulation time
sim_rel_tol = 1e-7;  % simulation accuracy.
sim_step = sim_time/1e4; % % simulation max step size.

% Load a network from the MatPower database:
switch choose_network
    case 1
        filename = 'case9'; % Matpower test-case network
        mpc = loadcase(filename); % Matpower function
    case 2
        filename = 'case_ieee30';
        mpc = loadcase(filename);
    case 3
        filename = 'case57';
        mpc = loadcase(filename);
        % remove duplicate branches in the network
        mpc.branch(20,:)=[];
        mpc.branch(35,:)=[];
    case 4
        % 118-bus network
        filename = 'case118';
        mpc = loadcase(filename);
        % remove duplicate branches
        mpc.branch(67,:)=[];
        mpc.branch(75,:)=[];
        mpc.branch(97,:)=[];
        mpc.branch(84,:)=[];
        mpc.branch(120,:)=[];
        mpc.branch(134,:)=[];
        mpc.branch(136,:)=[];
    case 5
        filename = 'case4gs';
        mpc = loadcase(filename);
    case 6
        filename = 'case6ww';
        mpc = loadcase(filename);
    case 7
        filename = 'case9Q';
        mpc = loadcase(filename);
    case 8
        filename = 'case30Q';
        mpc = loadcase(filename);
    otherwise
        disp('Unknown Matpower test-case');
        return
end

% Process and prepare network data:
% If the per-unit base voltage is not specified, choose a default value
ind = find(mpc.bus(:,BASE_KV)==0);
mpc.bus(ind,BASE_KV)=100; % [kV]
% Set a realistic per-unit base voltage for generators
mpc.bus(mpc.gen(:,GEN_BUS),BASE_KV)=12; % [kV]

% Make sure reference bus (slack bus) is bus no. 1
refbus = find(mpc.bus(:,BUS_TYPE) == 3); % reference bus
if (length(refbus)>1)
    disp('Error - more than one reference bus');
    return;
end
if (refbus~=1)
% Switch bus indices such that reference bus is bus no. 1
    mpc = switchBusInd(mpc,refbus,1);
end

% Sort bus data according to bus index
[~,ind1] = sort(mpc.bus(:,BUS_I));
mpc.bus = mpc.bus(ind1,:);
% Sort generator data according to generator bus index
[~,ind1] = sort(mpc.gen(:,GEN_BUS));
mpc.gen = mpc.gen(ind1,:);
% Define indices
ssbus = mpc.gen(:,GEN_BUS);  % generator bus indices (including ref bus)
ssgen = ssbus(2:end); % generator bus indices
% ssgen does NOT include the ref bus
% ssgen includes both the synchrnous machines and the photovoltaic inverter
if (synrp > length(ssgen))
    disp('Error - synrp variable exceeds the number of generators in the system'); beep; return;
end

% Shift voltage angles such that the angle at bus 1 is zero
mpc.bus(:,VA) = mpc.bus(:,VA) - mpc.bus(1,VA); 

% Per-unit definitions
baseVA = mpc.baseMVA * 1e6;  % [W] 1p.u. base power
baseV = 1e3*mpc.bus(:,BASE_KV); % [V] 1p.u. base voltage
baseA = baseVA./baseV; % [A] 1p.u. base current
baseZ = baseV./baseA; % [ohm] 1p.u. base impedance
ws = 2*pi*50; % system frequency [rad/sec]
N = size(mpc.bus,1);  % number of buses in the network

% Minimal values for network branches
% place minimal resistances on network branches (important for stability)
if ~exist('miniminal_branch_resistance_pu') 
    miniminal_branch_resistance_pu = 0.01;
end
mpc.branch(:,BR_R) = max( mpc.branch(:,BR_R), miniminal_branch_resistance_pu ); % branch resistance p.u.
% Place minimal charging susceptance on transmission lines 
% This is done to better represent the true dynamics of the transmission lines,
% and also to enable connection of syncronous machines in simulation
ind1 = find(mpc.branch(:,BR_B));
if isempty(ind1)
    zcharpu = 2.3; % [p.u.] characteristic impedance of transmission lines
else
    zcharpu = max((2*mpc.branch(ind1,BR_X)./mpc.branch(ind1,BR_B)).^0.5);
end
mpc.branch(:,BR_B) = max( mpc.branch(:,BR_B) , 2*mpc.branch(:,BR_X)/zcharpu^2 );

% Run power-flow - find the steady-state
mpc = runpf(mpc); 
clc;
if (mpc.success ~= 1)
    fprintf('%s\n\n',filename); % display network name
    disp('power flow failed to converge.');
    return
end

% Infinite bus parameters for simulation
vinf_d = (2^0.5)*mpc.bus(1,VM)*baseV(1); % [V] voltage source d-component
vinf_q = 0;
vinf_0 = 0;

% Synchronous machines basic parameters:
Prt = 2.5*mpc.gen(2:end,PG)*1e6; % [W] machine rated power (single phase) (maximum power for normal operation)
Prt = max(Prt,0.2*baseVA); % set minimal value for generators rated power
droop_Pref = mpc.gen(2:end,PG)*1e6; % [W] reference power for droop control (nominal generator power, single phase)
Ert = mpc.bus(ssgen,VM).*baseV(ssgen); % [Vrms] machine no-load output voltage at frequency w=ws
% Delete the machine being replaced with a photovoltaic inverter
Prt = Prt( [1:(synrp-1),(synrp+1):end]  );
Ert = Ert( [1:(synrp-1),(synrp+1):end]  );
droop_Pref = droop_Pref( [1:(synrp-1),(synrp+1):end]  );
% Other default parameters are based on this data:
poles = 2; % number of machine poles (must be even)
J = 0.005*Prt/ws; % rotor moment of inertia
poles_2Jws = (poles/2)^2./(J*ws); % constant used in simulation
Lpu = Ert.^2./(Prt*ws); % [H] per-unit inductance
Ld = 0.6*Lpu; % [H] direct axis synchronous inductance
Lq = Ld; % [H] quadrature axis synchronous inductance
L0 = 0.1*Ld; % [H] zero sequence inductance
Ls = (Ld+Lq)/2; % [H] synchronous inductance (used in simplified sync. machine model)
Ra = 0.05*Ert.^2./Prt; % [Ohm] armature resistance
Laf = (2^0.5)*5*Ert.^2./(ws*Prt); % [H] stator to rotor mutual inductance
Lff = 2*Laf.^2./Ld; % [H] field winding self inductance
Rf = Lff/1.2; % [Ohm] field winding resistance
If_dc = 0.2*Prt./Ert; % [A] field winding nominal DC current
Vf_dc = If_dc.*Rf; % [V] field winding nominal DC voltage
Lb2 = 2*Ld.*Lff - 3*Laf.^2; % constant used in simulation
% Droop control paramters
droop_D = 5e-4*(2/poles)./J; % [rad/sec / W] droop control sloop parameter

% % %%%% alternative version:
% poles = 2; % number of machine poles (must be even)
% J = 0.005*Prt/ws; % rotor moment of inertia
% poles_2Jws = (poles/2)^2./(J*ws); % constant used in simulation
% Vf_dc = 500 * ones(size(Ert)); % [V] field winding nominal DC voltage
% If_dc = Prt./(60*Vf_dc); % [A] field winding nominal DC current
% Rf = Vf_dc./If_dc; % [Ohm] field winding resistance
% Lpu = Ert.^2./(Prt*ws); % [H] per-unit inductance
% Ld = 0.3*Lpu; % [H] direct axis synchronous inductance
% Lq = Ld; % [H] quadrature axis synchronous inductance
% L0 = 0.1*Ld; % [H] zero sequence inductance
% Ls = (Ld+Lq)/2; % [H] synchronous inductance (may be used in simplified sync. machine model)
% Ra = 0.01*Ert.^2./Prt; % [ohm] armature resistance
% Laf = (2^0.5)*Ert./(ws*If_dc); % [H] stator to rotor mutual inductance
% Lff = 2*Laf.^2./Ld; % [H] field winding self inductance
% Lb2 = 2*Ld.*Lff - 3*Laf.^2; % constant used in simulation
% % Droop control paramters
% droop_D = 20./Prt; % [rad/sec / W] droop control sloop parameter

% Photovoltaic inverter parameters
% Rated values:
Pinv_rt = 0.1 *  3*mpc.gen(synrp+1,PG)*1e6; % [W] inverter rated power (total for three phase)
Pinv_rt  = max(Pinv_rt , 3e6); % [W] minimal value for rated power
Pdc_op = 0.6*Pinv_rt; % [W] inverter input power (total for three phase)
Vinv_rt = mpc.bus(ssgen(synrp),VM).*baseV(ssgen(synrp)); % [V] inverter rated output voltage (line-to-neutral RMS voltage)
% The following parameters are assigned default values that represent
% a possible working design. All parameters are calculated based on
% the rated power and voltage, and will scale properly for different rated values.
% These values can be used "as is" for preliminary system testings,
% or can be modified for a more accurate analysis.
Qstar = 0; % [VAr] desired reactive power (per-phase)
Linv = 10*(1/ws)*Pinv_rt/Vinv_rt^2; % [H] total inductance (per phase) at inverter output
Ebus = 0.4e-3*Pinv_rt; % [J] total bus capacitor stored energy at nominal bus voltage
duty_max = 0.9; % theoretical maximal duty-cycle required to generate the nominal line-to-ground voltage at each phase.
% current loop compensator constants:
Kcp = 10*4*pi*Linv; % proportional (P) constant
Kci = Kcp^2/(2*Linv); % integral (I) constant
% active power loop compensator constants:
Kpp = 200*(Ebus/Vinv_rt); % proportional (P) constant
Kpi = 0.1*Kpp^2/(Ebus/Vinv_rt); % integral (I) constant
% constants for inverter model:
Vdcset = (2^0.5)*Vinv_rt/duty_max; % [V] desired bus voltage (transformed to bridge input)
Cbuseff = 2*Ebus/Vdcset^2; % [F] effective bus capacitor (transformed to bridge input)

% Create a new database for constructing the transmission network dynamic model
mpcC = line2shunt(mpc); % converts line charging capacitors to shunt elements
% Calculate capacitors to be placed at the generators outputs in simulation:
Cbus_ext = 1e6*mpcC.bus(ssgen,BS)./(ws*baseV(ssgen).^2); % [F]
% mpcC.bus(ssbus,BS)=0; % remove shunt capacitors from generator buses
mpcC.bus(:,BS)=0; % remove all shunt capacitors from transmission network

% Extract load parameters from the MatPower database
% Loads are represented as series Rl shunt elements
% For details refer to function 'balloadsRL'
vv = baseV.*mpcC.bus(:,VM); % [V]
pp = 1e6*mpcC.bus(:,PD);  qq = 1e6*mpcC.bus(:,QD); % [VA]
ss2 = pp.^2 + qq.^2; % [VA]
if (any(qq<0))
    disp('Error: loads with negative reactive power are not supported');
    return;
end
if (any(pp<0))
    disp('Error: loads with negative active power are not supported');
    return;
end
R_bus = pp.*(vv.^2) ./ ss2;  % [Ohm] vector of load resistors
L_bus = (1/ws)*qq.*(vv.^2) ./ ss2; % [H] vector of load inductors
L_bus(ss2 == 0) = inf;   R_bus(ss2 == 0) = 0; % remove non-existing loads

% Build dynamic model of the network including loads using dq0 variables
% Build state-space model of the transmission network
[A1, B1, C1, D1, ~] = ssNetwMatPower(mpcC, ws);
% Build state-space model of loads:
[A2, B2, C2, D2] = balLoadsRL(R_bus, L_bus, ws);
% Merge state-space models of original network and load network:
[A3, B3, C3, D3] = mergeParlNetw(A1, B1, C1, D1, A2, B2, C2, D2);
% Eliminate all buses except from generator buses:
[Anet, Bnet, Cnet, Dnet] = elimBus(A3, B3, C3, D3, ssbus);

clear A1 B1 C1 D1;
clear A2 B2 C2 D2;
clear A3 B3 C3 D3;

% Compute initial capacitor voltages for simulation - using the power flow solution
vppa=(2^0.5)*baseV(ssgen).*exp(1j*mpcC.bus(ssgen,VA)*pi/180);
vcap_init_d = real(vppa); % [V]
vcap_init_q = imag(vppa); % [V]
vcap_init_0 = 0*vppa; % [V]

% Signal routing for simulation
nG = length(ssgen); % number of generators not including bus 1
units_to_net = [1,4:(3+nG),2,(4+nG):(3+2*nG),3,(4+2*nG):(3+3*nG)];
[~,net_to_units] = sort(units_to_net);
routing1 = [1:(synrp-1),(synrp+1):nG];
routing1 = [routing1 , routing1+nG, routing1+2*nG];
routing2 = [synrp , synrp+nG , synrp+2*nG];
temp = [routing1 , routing2];
[~,routing3] = sort(temp);

% A short status report
fprintf('%s\n',filename);
fprintf('Number of buses = %d\n\n',N);
fprintf('Please wait ...\n');
fprintf('Simulation time = %d sec\n\n',round(sim_time));
if (N>100)
    fprintf('Simulation running...\n');
    fprintf('Warning: long simulation time. This may take several minutes...\n');
elseif (N>20)
    fprintf('Simulation running...\n');
end

% Run simulation (Simulink)
fprintf('Starting Simulink...\n\n');
exSyncGenPVinv_sim;
set_param(bdroot,'StopTime',num2str(sim_time));
set_param(bdroot,'RelTol',num2str(sim_rel_tol));
set_param(bdroot,'MaxStep',num2str(sim_step));
sim(bdroot);

%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Visualize results:
%%%%%%%%%%%%%%%%%%%%%%%%
set(0,'defaulttextinterpreter','latex')
set(0,'defaultfigurecolor',[1 1 1])
set(0,'defaultaxesfontsize',9);

figure(1);
numplots = 1;
cur_plot = 1;
xmin = 1; % start time for plot
xmax = sim_time; % final time for plot

subplot(numplots,1,cur_plot);
plot(tsim,netp,'-k','Color',[0 0 0],'LineWidth',0.7);
xlim([xmin xmax]);
ylabel('$P$ [MW]','FontSize',9);
cur_plot = cur_plot+1;

xlabel('Time [s]');

axesHandles = findall(0,'type','axes');
set(axesHandles,'TickLabelInterpreter', 'latex')
