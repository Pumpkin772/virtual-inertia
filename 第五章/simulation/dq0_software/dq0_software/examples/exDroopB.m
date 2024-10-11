% Dynamics and stability analysis of large power systems
% with physical synchronous machines.
% Short description:
% a) Droop control
% b) No infinite bus (steady-state frequency determined by control)
% c) DQ0 models based on a central angle reference frame
%
% The user may choose between transient simulation
% or small signal analysis.
%
% ** Associated Simulink file: 'exDroopB_sim'
%
% To begin, choose an example below, and run the script.
% Note: the MatPower package should be installed (available online).

clc;
close all;

% Choose a test-case system from the list below:
choose_system = 1;

% Choose task:
% 1 = transient simulation using the (physical) non-linear model
% 2 = linear analysis. Linearize the model at the operating point,
% and compute poles and frequency response
choose_task = 1;

% Simulation parameters
sim_time = 15;  % [s] simulation time
sim_rel_tol = 1e-8;  % simulation accuracy
sim_step = sim_time/1e4; % simulation max step size
% Define input steps for transient analysis:
% (ignored in small-signal analysis)
step_gen_num = 1; % choose generator to apply step
step_time = 10; % [s]
step_size = 20e6; % [W]

define_constants; % MatPower function. Defines MatPower constants.

% Load a network from the MatPower database:
switch choose_system
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
mpc.bus(mpc.gen(:,GEN_BUS),BASE_KV)=20; % [kV]

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
ssgen = ssbus; % generator buses. In this example - including ref bus.

% Shift voltage angles such that the angle at bus 1 is zero
mpc.bus(:,VA) = mpc.bus(:,VA) - mpc.bus(1,VA); 

% Per-unit definitions
baseVA = mpc.baseMVA * 1e6;  % [W] 1p.u. base power
baseV = 1e3*mpc.bus(:,BASE_KV); % [V] 1p.u. base voltage
baseA = baseVA./baseV; % [A] 1p.u. base current
baseZ = baseV./baseA; % [Ohm] 1p.u. base impedance
ws = 2*pi*50; % system frequency [rad/s]
N = size(mpc.bus,1);  % number of buses in the network

% Minimal values for network branches
% Place minimal resistances on network branches (important for stability)
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

% Run power-flow
mpc = runpf(mpc); 
if (mpc.success ~= 1)
    fprintf('%s\n\n',filename); % display network name
    disp('power flow failed to converge.');
    return
end

% Synchronous machines parameters
% Basic parameters
Prt = 2.5*mpc.gen(1:end,PG)*1e6; % [W] machine rated power (single phase) (maximum power for normal operation)
Prt = max(Prt,0.2*baseVA); % set minimal value for generators rated power
Ert = mpc.bus(ssgen,VM).*baseV(ssgen); % [Vrms] machine no-load output voltage at frequency w=ws
J_GAIN_CNST = 0.005; % [s^2] constant defining moment of inertia
KD_GAIN_CNST = 0; % 100; % constant defining damping factor
DROOP_R_GAIN_CNST = 10; % [rad/s] constant defining droop control sloop parameter
% Other default parameters are based on this data:
Pm_rt = 3*Prt; % [W] rated mechanical (input) power
poles = 2; % number of machine poles (must be even)
J = J_GAIN_CNST*Prt/ws; % rotor moment of inertia
Jtot = sum(J);
Kd = KD_GAIN_CNST*J*ws; % damping factor
poles_2Jws = poles./(2*J*ws); % constant used in simulation
Vf_dc = 500 * ones(size(Ert)); % [V] field winding nominal DC voltage
If_dc = Prt./(60*Vf_dc); % [A] field winding nominal DC current
Rf = Vf_dc./If_dc; % [Ohm] field winding resistance
Lpu = Ert.^2./(Prt*ws); % [H] per-unit inductance
Ld = 0.3*Lpu; % [H] direct axis synchronous inductance
Lq = Ld; % [H] quadrature axis synchronous inductance
L0 = 0.1*Ld; % [H] zero sequence inductance
Ra = 0.01*Ert.^2./Prt; % [ohm] armature resistance
Laf = (2^0.5)*Ert./(ws*If_dc); % [H] stator to rotor mutual inductance
Lff = 2*Laf.^2./Ld; % [H] field winding self inductance
Lb2 = 2*Ld.*Lff - 3*Laf.^2; % constant used in simulation
Ls = (Ld+Lq)/2; % [H] synchronous inductance (may be used in simplified sync. machine model)
% Droop control paramters
droop_Pref = mpc.gen(:,PG)*1e6; % [W] reference power for droop control (nominal generator power, single phase)
droop_R = DROOP_R_GAIN_CNST./Prt; % [rad/s/W] droop control sloop parameter

% Create a new database for constructing the transmission network dynamic model
mpcC = line2shunt(mpc); % converts line charging capacitors to shunt elements
% Calculate capacitors to be placed at the sync machine outputs in simulation:
Cbus_ext = 1e6*mpcC.bus(ssgen,BS)./(ws*baseV(ssgen).^2); % [F]
% mpcC.bus(ssbus,BS)=0; % remove shunt capacitors from generator buses
mpcC.bus(:,BS)=0; % remove all shunt capacitors from transmission network

% Extract load parameters from the MatPower database
% Loads are represented as series RL shunt elements.
% For details refer to function 'balLoadsRL'
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
vppa=(2^0.5)*mpc.bus(ssgen,VM).*baseV(ssgen).*exp(1j*mpc.bus(ssgen,VA)*pi/180);
vcap_init_d = real(vppa); % [V]
vcap_init_q = imag(vppa); % [V]
vcap_init_0 = 0*vppa; % [V]

% Short status report
fprintf('%s\n',filename);
fprintf('Number of buses = %d\n\n',N);
fprintf('Simulation time = %d s\n\n',round(sim_time));

% Run simulation (Simulink)
fprintf('Starting Simulink...\n\n');
exDroopB_sim;
set_param(bdroot,'StopTime',num2str(sim_time));
set_param(bdroot,'RelTol',num2str(sim_rel_tol));
set_param(bdroot,'MaxStep',num2str(sim_step));

% Small-signal analysis
if (choose_task == 2)
    fprintf('Small signal analysis... Please wait ...\n');
    if (N>100)
        fprintf('Simulation converging to steady state.\n');
        fprintf('Warning: long simulation time. This may take several minutes...\n');
    elseif (N>20)
        fprintf('Simulation converging to steady state...\n');
    end
% Cancel input steps:
    step_size = 0;
% Simulate model and linearize at operating point:
    io = getlinio(bdroot);
    tic
    linsys = linearize(bdroot,io,sim_time);
    toc
    [Ap,Bp,Cp,Dp]=ssdata(linsys); % linear system
% The linear model inputs and outputs are marked in Simulink.
% Input pertubation at droop-control referenece power
% output measuement at transmission network dq0 currents
    eeg = eig(Ap); % compute eigenvalues
% State matrix of linear system
    figure(1);
    spy(Ap);
    title('state matrix A of linear system');
% Pole map
    figure(2);
    plot(real(eeg),imag(eeg),'kx');
    xlim([-100 10]);
    xlabel('Real Axis');
    ylabel('Imaginary Axis');
    title('Pole map');
% Bode plots
    sys_inputs = 1;
    sys_outputs = 1:2;
    figure(3);
    bode(linsys(sys_outputs,sys_inputs),{1e-2,1e5});
    return
end

% Transient simulation
fprintf('Transient simulation... Please wait ...\n')
if (N>100)
    fprintf('Warning: long simulation time. This may take several minutes...\n');
elseif (N>20)
    fprintf('Simulation converging to steady state...\n');
end
tic
sim(bdroot);
toc

%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Visualize results:
%%%%%%%%%%%%%%%%%%%%%%%%
% Time-domain
set(0,'defaulttextinterpreter','latex')
set(0,'defaultfigurecolor',[1 1 1])
set(0,'defaultaxesfontsize',9);
figure(1);
numplots = 3;
cur_plot = 1;
xmin = step_time-0.5; % start time for plot
xmax = step_time+1; % final time for plot

subplot(numplots,1,cur_plot);
plot(tsim,p/1e6,'-k','Color',[0.25,0.25,0.25],'LineWidth',0.7);
xlim([xmin xmax]);
ylabel('$P$ [MW]','FontSize',9);
hold off; cur_plot = cur_plot+1;

subplot(numplots,1,cur_plot);
plot(tsim,f,'-k','Color',[0.25,0.25,0.25],'LineWidth',0.7);
hold on;
plot(tsim,wc/(2*pi),'--k','Color',[0.25,0.25,0.25],'LineWidth',0.7);
xlim([xmin xmax]);
ylabel('$f$ [Hz]','FontSize',9);
hold off; cur_plot = cur_plot+1;

subplot(numplots,1,cur_plot);
plot(tsim,delta,'-k','Color',[0.25,0.25,0.25],'LineWidth',0.7);
xlim([xmin xmax]);
ylabel('$\delta$ [deg]','FontSize',9);
hold off; cur_plot = cur_plot+1;

xlabel('Time [s]');

axesHandles = findall(0,'type','axes');
set(axesHandles,'TickLabelInterpreter', 'latex')
