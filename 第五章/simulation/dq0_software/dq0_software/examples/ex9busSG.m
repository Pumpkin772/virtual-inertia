% Simulate the dynamics of a 9-bus power system in the time domain,
% using dq0 signals. The simulation includes dynamic
% models of the network, loads, and synchronous generators.

% ** To start the simulation - run this script.
% ** Associated Simulink file: 'ex9busSG_sim'

% The system in this example includes three major blocks:
% the transmission network, the loads, and the generators.
% The complete system is simulated in the
% time domain, using Simulink.

% The system converges to a steady-state that matches
% its power flow solution. The user may trigger a dynamic
% response by defining an input step in mechanical power.

% For more information about dq0 signals and modeling,
% refer to the documentation in function 'ssNetw'

clc;
close all;

% Load data of 9-bus test-case network
load('ex9busSG_data');

% Simulation parameters
sim_time = 6;  % [s] simulation time
sim_rel_tol = 1e-5;  % simulation accuracy
sim_step = sim_time/1e3; % simulation max step size

% Define step in mechanical power input:
which_generator_to_step = 2; % index of gen.
ppm_step_time = 2; % [s] when the step occurs
ppm_rel_step = 0.1; % step size relative to generator max. power
% To cancel the step define: ppm_rel_step=0
% In this case the system will converge to an operating point
% predicted by its power-flow solution

ws = 2*pi*50; % [rad/s] nominal system frequency
N = size(mpc.bus,1); % number of buses in the network

% Per-Unit quantities
baseVA = mpc.baseMVA*1e6;  % [W] 1p.u. base power
baseA = baseVA/baseV; % [A] 1p.u. base current
baseZ = baseV/baseA; % [Ohm] 1p.u. base impedance

% Buses of generators:
Genbus = mpc.gen(:,GEN_BUS);  % generator bus indices (not including bus 1)
if (Genbus(1)==1)
    Genbus=Genbus(2:end);
end % remove bus 1
ssbus = [1 ; Genbus];  % bus indices of generators (including bus 1)

% Build state-space model of the network and loads
% Build state-space model of original network:
[A1, B1, C1, D1, YbusPU] = ssNetwMatPower(mpc, ws); 
% Extract load parameters from the database:
% Loads are represented as series RL shunt elements
% For details refer to function 'balLoadsRL'
vv = 1e3*mpc.bus(:,BASE_KV).*mpc.bus(:,VM);
pp = 1e6*mpc.bus(:,PD);   qq = 1e6*mpc.bus(:,QD);
ss2 = pp.^2 + qq.^2;
if (any(qq<0))
    disp('Warning: negative load reactive power is not supported');
end
if (any(pp<0))
    disp('Warning: negative load active power is not supported');
end
R_bus = pp.*(vv.^2) ./ ss2;  % vector of load resistors
L_bus = (1/ws)*qq.*(vv.^2) ./ ss2; % vector of load inductors
L_bus(ss2 == 0) = inf;   R_bus(ss2 == 0) = 0; % remove non-existing loads
% Build state-space model of loads:
[A2, B2, C2, D2] = balLoadsRL(R_bus, L_bus, ws);
% Merge state-space models of original network and load network:
[A, B, C, D] = mergeParlNetw(A1, B1, C1, D1, A2, B2, C2, D2);
% Eliminate all buses except from generator buses:
[Ar, Br, Cr, Dr] = elimBus(A, B, C, D, ssbus);
Ar=full(Ar); Br=full(Br);
Cr=full(Cr); Dr=full(Dr); % convert matrices to full for Simulink

% Simulink  - parameters for infinite bus (bus 1)
VV1ph = (mpc.bus(1,VM)*baseV) * exp(1i*mpc.bus(1,VA)*pi/180);
ref_vd_MKS = (2^0.5)*real(VV1ph);
ref_vq_MKS = (2^0.5)*imag(VV1ph);

% Simulink -  parameters of synchronous machines
Pn_gens = mpc.gen((1:length(Genbus))+1,PG)*1e6;  % nominal output power one phase
Qn_gens = mpc.gen((1:length(Genbus))+1,QG)*1e6;  % nominal output power one phase
Vn_gens = mpc.bus(Genbus,VM)*baseV; % nominal voltage magnitude of one phase
delta_init = mpc.bus(Genbus,VA)*pi/180; % initial power angle
ppm_default = (3*Pn_gens);  % mechanical power input 
Pmax_gens = 12*Pn_gens;  % max. generator power output of three-phase
Pmax_gens = max(Pmax_gens, 0.1*sum(mpc.gen(:,PG))*1e6 );
ppm_step_size_W = ppm_rel_step*Pmax_gens(which_generator_to_step);
poles = 2;
Ls = 3*(Vn_gens.^2) ./ (ws*Pmax_gens);
lamda = (2^0.5) * Vn_gens/ws;
J_GAIN_CNST = 0.002; % default=0.002
J = J_GAIN_CNST*Pmax_gens/ws;
KD_GAIN_CNST = 100;% [1/s] default=100
Kd = KD_GAIN_CNST*J*ws;
poles_2Jws = poles./(2*ws*J);
wsLs = ws*Ls;

% Simulink - signal routing data
MM = length(ssbus);
gvec = 1+(1:length(Genbus))';
Bus_to_Ref = [1; 1+MM; 1+2*MM];
Bus_to_Gen = [gvec; gvec+MM; gvec+2*MM];
[temp,ind] = sort([Bus_to_Ref; Bus_to_Gen]);
Units_to_Net = ind;
% Signal routing for mechanical power input step:
if (which_generator_to_step>length(ppm_default))
    disp('Error: invalid value for ''which_generator_to_step''');
    return
end
[temp,ind] = sort([which_generator_to_step, ...
    1:(length(ppm_default)-1)] );
step_index_vector = ind;
 
% Run simulink
ex9busSG_sim;
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
tmin = 1;
tmax = sim_time;

subplot(3,1,1);
plot(tsim, bus1_current_abc/baseA,'LineWidth',0.5);
xlim([tmin tmax]);
ylabel('$I_{abc,1}$ [p.u.]');

subplot(3,1,2);
plot(tsim, bus1_current_dq0(:,1)/baseA,'-k','Color',[0 0 0],'LineWidth',0.6);
xlim([tmin tmax]);
ylabel('$I_{d,1}$ [p.u.]');

subplot(3,1,3);
plot(tsim, bus1_current_dq0(:,2)/baseA,'-k','Color',[0 0 0],'LineWidth',0.6);
xlim([tmin tmax]);
ylabel('$I_{q,1}$ [p.u.]');

xlabel('Time [s]');

axesHandles = findall(0,'type','axes');
set(axesHandles,'TickLabelInterpreter', 'latex')
