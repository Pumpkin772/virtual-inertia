% Dynamic simulation of the IEEE 14 bus test case network, with
% physical synchronous machine models and photovoltaic inverters (renewable sources).
% Simulation is based on dq0 quantities.
% The user may choose between transient simulation or small signal analysis.
%
% Note: to run this script, the MatPower package
% should be installed (available online).
%
% ** Associated Simulink file: 'ex14busSGandPVSim'
%
% Network structure:
% the 14-bus network topology and parameters are
% loaded from the Matpower database.
% bus 1 - an infinite bus (voltage source)
% bus 2 - synchronous generator
% bus 3 - synchronous generator
% bus 6 - photovoltaic inverter
% bus 8 - photovoltaic inverter
%
% The synchronous machine model is based on a
% book "Electric Machinery" by Fitzgerald.
% We represent input voltages and output currents in
% a references frame rotating with angle ws*t,
% (instead of theta - the rotor electrical angle),
% where ws is the frequency of the infinite bus.
% This allows direct connection of the machine model
% to the network.

% The models of the photovoltaic generators (inverters)
% are based on a paper "Bus Voltage Control With Zero
% Distortion and High Bandwidth for Single-Phase Solar Inverters"

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear all;

% Please choose a task:
% 1 = transient simulation using the (physical) nonlinear model
% 2 = small signal analysis. Compute operating point,
% derive a linear small signal model, and compute
% state-matrix eigenvalues
choose_task = 1;

% Simulation parameters
sim_time = 8;  % [s] simulation time
sim_rel_tol = 1e-4;  % simulation accuracy
sim_step = sim_time/1e4; % % simulation max step size
% Define input steps for transient analysis:
% (ignored in small-signal analysis)
step_time_bus2 = 3; % [s]
step_size_bus2 = 3*20e6; % [W] typical value 3*20e6
% mechnical power synchronous generator bus 3
step_time_bus3 = 3; % [s]
step_size_bus3 = 0;  % [W] typical value 3*10e6
% DC input power photovoltaic inverter bus 6
step_time_bus6 = 5; % [s]
step_size_bus6 = 3*5e6; % [W] typical value 3*5e6
% DC input power photovoltaic inverter bus 8
step_time_bus8 = 5; % [s]
step_size_bus8 = 0; % [W] typical value 3*5e6

% Load network Matpower database
filename = 'case14'; % IEEE 14-bus test-case network
mpcA = loadcase(filename);

% Prepare network data:
define_constants; % MatPower function. Defines MatPower constants.
ws = 2*pi*50; % [rad/s] nominal system frequency (infinite bus frequency)
NN = size(mpcA.bus,1);  % number of buses in the network
mpcA = rmfield(mpcA,'gencost'); mpcA=rmfield(mpcA,'bus_name'); % remove unused fields
% Per-unit base voltage default value
ind = find(mpcA.bus(:,BASE_KV)==0);
mpcA.bus(ind,BASE_KV)=100; % [kV]
% Place minimal resistances on network branches (important for stability)
if ~exist('miniminal_branch_resistance_pu') 
    miniminal_branch_resistance_pu = 0.01;
end
mpcA.branch(:,BR_R) = max( mpcA.branch(:,BR_R), miniminal_branch_resistance_pu ); % branch resistance p.u.
% Eliminate loads with negative reactive power
mpcA.bus(  mpcA.bus(:,QD) < 0  ,QD) = 0;
% Remove all shunt capacitors except on buses 2,3
mpcA = line2shunt(mpcA);
mpcA.bus(1,BS)=0; mpcA.bus(4,BS)=0;
mpcA.bus(5,BS)=0; mpcA.bus(9,BS)=0;
% Run power flow
mpcA = runpf(mpcA);  clc;
% Define per-unit (p.u.) quantities and other constants
baseVA = mpcA.baseMVA * 1e6;  % [W] 1p.u. base power
baseV = 1e3*mpcA.bus(1,BASE_KV); % [Vrms] 1p.u. base voltage
baseA = baseVA/baseV; % [A] 1p.u. base current
baseZ = baseV/baseA; % [ohm] 1p.u. base impedance

% Synchronous machine(s) parameters
% (relates to synchronous generators on buses 2,3)
% enter the following specifications:
Prt = [100 ; 50]*1e6; % [W] machine rated power (single phase) (maximum power for normal operation)
Pgm = 3*[40; 20]*1e6; % [W] machine mechanical power input at steady-state
% (equals to power output of three phases not including loses)
Ert = mpcA.gen([2 3],VG)*baseV; % [Vrms] machine no-load output voltage at frequency w=ws
% Default parameters based on this data 
Pm_rt = 3*Prt; % [W] rated mechanical (input) power
poles = 2; % number of machine poles (must be even)
J_GAIN_CNST = 0.005; % [s^2] constant defining moment of inertia
KD_GAIN_CNST = 10; % constant defining damping factor
J = J_GAIN_CNST*Prt/ws; % rotor moment of inertia
Kd = KD_GAIN_CNST*J*ws; % damping factor
poles_2Jws = poles./(2*J*ws); % constant used in simulation
If_dc = 0.2*Prt./Ert; % [A] field winding nominal DC current
Lpu = Ert.^2./(Prt*ws); % [H] per-unit inductance
Ld = 0.6*Lpu; % [H] direct axis synchronous inductance
Lq = Ld; % [H] quadrature axis synchronous inductance
L0 = 0.1*Ld; % [H] zero sequence inductance
Ra = 0.05*Ert.^2./Prt; % [ohm] armature resistance
Laf = (2^0.5)*Ert./(ws*If_dc); % [H] stator to rotor mutual inductance
Lff = 2*Laf.^2./Ld; % [H] field winding self inductance
Rf = Lff/1.2; % [ohm] field winding resistance
Vf_dc = If_dc.*Rf; % [V] field winding nominal DC voltage
Lb2 = 2*Ld.*Lff - 3*Laf.^2; % constant used in simulation
Ls = (Ld+Lq)/2; % [H] synchronous inductance (may be used in simplified sync. machine model)

% Photovoltaic generator parameters
% (relates to generators on buses 6,8)
Pren = [10 ; 10]*1e6; % [W] steady-state output power of photovoltaic source
% for a single phase (total power is 3 times higher)
Cbus = [600e-6 ; 600e-6]; % [F] bus capacitor (at the inverter input)
Vbus_ref = [800 ; 800]; % [V] bus voltage reference
kp = [1.3e-6 ; 1.3e-6];  % proportional constant of inverter PI controller 
ki = [0.3 ; 0.3]; % integral constant of inverter PI controller
Cinv_pu = [0.1 ; 0.1]; % [p.u.] photovoltaic inverter output capacitance, units:
% susceptance in p.u. (equals ws*C in p.u.)
Cinv =  Cinv_pu/(ws*baseZ); % [F] inverter output capacitance

% The following code defines two network databases:
% mpcC - is used for constructing a dynamic model of the transmission network
% mpcE - is used for power flow calculations
% Prepare network database for constructing a dynamic model
mpcC = mpcA; % new database for constructing the transmission network dynamic model
mpcC = line2shunt(mpcC); % converts line charging capacitors to shunt elements
% Calculate capacitors to be placed at the sync machine output in simulation:
Cbus2 = (1e6*mpcC.bus(2,BS))/(ws*baseV^2);% [F] shunt capacitor on bus 2
Cbus3 = (1e6*mpcC.bus(3,BS))/(ws*baseV^2);% [F] shunt capacitor on bus 3
% Remove these capacitors from the transmission network database:
mpcC.bus(2,BS)=0;
mpcC.bus(3,BS)=0;

% Extract load parameters from the MatPower database
% Loads are represented as series L-R shunt elements
% For details refer to function 'balLoadsRL'
vv = baseV*mpcC.bus(:,VM); %[V]
pp = 1e6*mpcC.bus(:,PD);  qq = 1e6*mpcC.bus(:,QD); % [VA]
ss2 = pp.^2 + qq.^2; % [VA]
if (any(qq<0))
    disp('Warning: loads with negative reactive power are not supported');
    return;
end
if (any(pp<0))
    disp('Warning: loads with negative active power are not supported');
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
ssbus = [1; 2; 3; 6; 8]; % buses to keep
[Aksi, Bksi, Cksi, Dksi] = elimBus(A3, B3, C3, D3, ssbus);
Anf = full(Aksi); Bnf = full(Bksi);
Cnf = full(Cksi); Dnf = full(Dksi);

% Build routing vectors for Simulink
Nss = length(ssbus);
units_to_net = zeros(1,3*Nss);
for ii=1:Nss
    units_to_net(ii:Nss:end) = (ii-1)*3 + [1 2 3];
end
[~,net_to_units] = sort(units_to_net);

% Prepare network database for power flow calculations
mpcE = mpcC; % new database for power flow calculations
% Add two new buses 15,16 to represent the synchronous machines behind
% their synchronous impedances.
% Copy original data
mpcE.bus = zeros(16,size(mpcC.bus,2));
mpcE.bus(1:size(mpcC.bus,1),1:size(mpcC.bus,2)) = mpcC.bus;
% Define new buses
mpcE.bus([15; 16],BUS_I)= [15; 16]; % new bus indices
mpcE.bus([15; 16],BUS_TYPE)= [2; 2]; % define new buses as PV buses
mpcE.bus([2; 3],BUS_TYPE)= [1; 1]; % define original buses as PQ buses
mpcE.bus([15; 16],PD)= [0; 0]; % real power demand
mpcE.bus([15; 16],QD)= [0; 0]; % reactive power demand
mpcE.bus([15; 16],BUS_AREA)= [1; 1]; 
mpcE.bus([15; 16],VM) = Ert/baseV;
mpcE.bus([15; 16],BASE_KV)= baseV/1e3; 
mpcE.bus([15; 16],ZONE)= 1; 
mpcE.bus([15; 16],VMAX)= 1.06; 
mpcE.bus([15; 16],VMIN)= 0.94;
% Add two new branches to connect the new buses of the syncronous machines
% Copy original data
mpcE.branch = zeros(22,size(mpcC.branch,2));
mpcE.branch(1:size(mpcC.branch,1),1:size(mpcC.branch,2)) = mpcC.branch;
% Define new branches
mpcE.branch([21; 22],F_BUS) = [15; 16];
mpcE.branch([21; 22],T_BUS) = [2; 3];
mpcE.branch([21; 22],BR_R) = Ra/baseZ;
mpcE.branch([21; 22],BR_X) = (Ls*ws)/baseZ;
mpcE.branch([21; 22],BR_STATUS) = 1; % active line
mpcE.branch([21; 22],ANGMIN) = -360;
mpcE.branch([21; 22],ANGMAX) = +360;
% Define synchronous generators on the new buses
mpcE.gen = zeros(3,size(mpcC.gen,2)); 
mpcE.gen(1,:) = mpcC.gen(1,:); % copy infinite bus generator
mpcE.gen([2; 3],GEN_BUS) = [15; 16]; % new generators on buses 15,16
mpcE.gen([2; 3],PG) = Pgm/(3e6);
mpcE.gen(:,QMAX) = 100;
mpcE.gen(:,QMIN) = -100;
mpcE.gen([2; 3],VG) = Ert/baseV;
mpcE.gen([2; 3],MBASE) = 100;
mpcE.gen([2; 3],GEN_STATUS) = 1;
mpcE.gen([2; 3],PMAX) = 300;
mpcE.gen([2; 3],PMIN) = 0;
% Add photovoltaic generators as PQ sources on buses 6 and 8
mpcE.bus([6; 8],BUS_TYPE)= [1; 1]; % define buses as PQ buses
% Add shunt capacitors
mpcE.bus(2,BS) = mpcE.bus(2,BS) + Cbus2*ws*baseV^2/1e6;
mpcE.bus(3,BS) = mpcE.bus(3,BS) + Cbus3*ws*baseV^2/1e6;
mpcE.bus([6, 8],BS) = mpcE.bus([6, 8],BS) + (baseV^2)*ws*Cinv/1e6;
% Represent loads as shunt elements
R_bus16 = [R_bus; inf; inf];
L_bus16 = [L_bus; 0; 0];
mpcE.bus(:,GS) = mpcE.bus(:,GS) + (baseV^2)*real(1./(R_bus16-1j*ws*L_bus16))/1e6;
mpcE.bus(:,BS) = mpcE.bus(:,BS) - (baseV^2)*imag(1./(R_bus16-1j*ws*L_bus16))/1e6;
mpcE.bus(:,PD) = 0;
mpcE.bus([6, 8],PD) = -Pren/1e6; % add power generation of photovoltaic inverters as negative loads
mpcE.bus(:,QD) = 0;
% Run power flow
mpcE = runpf(mpcE);
if (mpcE.success == 0)
    disp('Power flow did not converge');
    return;
end

% Infinite bus parameters for simulation
vinf_d = (2^0.5)*mpcA.bus(1,VM)*baseV ; % [V] voltage source d-component
vinf_q = 0;
vinf_0 = 0;

% Compute initial bus voltages (for simulation) using the power flow solution
busii = 2; vpp=(2^0.5)*mpcE.bus(busii,VM)*baseV*exp(1j*mpcE.bus(busii,VA)*pi/180);
vinit_dq0_bus2 = [real(vpp) ; imag(vpp) ; 0]; % [Vd; Vq; V0], units: V
busii = 3; vpp=(2^0.5)*mpcE.bus(busii,VM)*baseV*exp(1j*mpcE.bus(busii,VA)*pi/180);
vinit_dq0_bus3 = [real(vpp) ; imag(vpp) ; 0]; % [Vd; Vq; V0], units: V
busii = 6; vpp=(2^0.5)*mpcE.bus(busii,VM)*baseV*exp(1j*mpcE.bus(busii,VA)*pi/180);
vinit_dq0_bus6 = [real(vpp) ; imag(vpp) ; 0]; % [Vd; Vq; V0], units: V
busii = 8; vpp=(2^0.5)*mpcE.bus(busii,VM)*baseV*exp(1j*mpcE.bus(busii,VA)*pi/180);
vinit_dq0_bus8 = [real(vpp) ; imag(vpp) ; 0]; % [Vd; Vq; V0], units: V


% Run simulation (Simulink)
disp('Starting Simulink...')
ex14busSGandPV_sim;
set_param(bdroot,'StopTime',num2str(sim_time));
set_param(bdroot,'RelTol',num2str(sim_rel_tol));
set_param(bdroot,'MaxStep',num2str(sim_step));

% Small-signal analysis
% Linearize the model at operating point and compute eigenvalues
if (choose_task == 2)
    disp('Small signal analysis... Please wait ...')
% Cancel input steps:
    step_size_bus2 = 0;
    step_size_bus3 = 0;
    step_size_bus6 = 0;
    step_size_bus8 = 0;
% Simulate model and linearize at operating point:
    io = getlinio(bdroot);
    linsys = linearize(bdroot,io,3.5);
    [Ap,Bp,Cp,Dp]=ssdata(linsys); % linear system
% The linear model inputs and outputs are marked in Simulink.
% The order is:
% input = [pm2; pm3; pdc6; pdc8]
% output = [delta2; p2; delta3; p3; vbus_6; pinv_6; vbus_8; pinv_8]
    eeg = eig(Ap); % compute eigenvalues
    fprintf('eEigenvalue maximal real part is %f\n',max(real(eeg)));
% Bode plots
    figure(1);
    bode(linsys);    
% Pole map
    figure(2);
    pzmap(linsys);
    xlim([-500 0]);
    return;
end

% Transient simulation
disp('Transient simulation... Please wait ...')
sim(bdroot);


%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Visualize results:
%%%%%%%%%%%%%%%%%%%%%%%%
set(0,'defaulttextinterpreter','latex')
set(0,'defaultfigurecolor',[1 1 1])
set(0,'defaultaxesfontsize',9);
figure(1);
numplots = 5;
cur_plot = 1;
xmin = 2; % start time for plot
xmax = sim_time; % final time for plot

subplot(numplots,1,cur_plot);
plot(tsim,pm2/1e6,'-k','Color',[0 0 0],'LineWidth',0.7);
xlim([xmin xmax]);
ylabel('$P_{m,2}$ [MW]','FontSize',9);
cur_plot = cur_plot+1;

subplot(numplots,1,cur_plot);
plot(tsim,p2/1e6,'-k','Color',[0 0 0],'LineWidth',0.7);
xlim([xmin xmax]);
ylabel('$P_2$ [MW]','FontSize',9);
cur_plot = cur_plot+1;

subplot(numplots,1,cur_plot);
plot(tsim,delta2,'-k','Color',[0 0 0],'LineWidth',0.7);
xlim([xmin xmax]);
ylabel('$\delta_2$ [deg]','FontSize',9);
cur_plot = cur_plot+1;

subplot(numplots,1,cur_plot);
plot(tsim,pinv_6/1e6,'-k','Color',[0 0 0],'LineWidth',0.7);
xlim([xmin xmax]);
ylabel('$P_{inv,6}$ [MW]','FontSize',9);
cur_plot = cur_plot+1;

subplot(numplots,1,cur_plot);
plot(tsim,ph6,'-k','Color',[0 0 0],'LineWidth',0.7);
xlim([xmin xmax]);
ylabel('$\theta_6$ [deg]','FontSize',9);
cur_plot = cur_plot+1;

xlabel('Time [s]');

axesHandles = findall(0,'type','axes');
set(axesHandles,'TickLabelInterpreter', 'latex')
