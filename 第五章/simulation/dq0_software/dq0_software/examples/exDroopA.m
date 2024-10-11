% Dynamics and stability analysis of large power systems
% with physical synchronous machines.
% Short description:
% a) Droop control
% b) Infinite bus
% c) DQ0 models based on the unified reference frame
%
% Demonstrates convegence to steady-state,
% and linearization at operating point.
%
% ** Associated Simulink file: 'exDroopA_sim'
%
% To begin, choose an example below, and run the script.
% Note: the MatPower package should be installed.
% (available online).

clc;
close all;

% Choose an example from the list below:
choose_example = 1;

define_constants; % MatPower function. Defines MatPower constants.

%%%%%%%%%% Simulation parameters %%%%%%%%%%%%%
sim_time = 30;  % [sec] simulation time
sim_rel_tol = 1e-7;  % simulation accuracy.
sim_step = sim_time/1e4; % % simulation max step size.

% Load a network from the MatPower database:
switch choose_example
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
ssgen = ssbus(2:end); % generator bus indices NOT includeing ref bus

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

% Run power-flow
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

% Synchronous machines parameters basic parameters:
Prt = 2.5*mpc.gen(2:end,PG)*1e6; % [W] machine rated power (single phase) (maximum power for normal operation)
Prt = max(Prt,0.2*baseVA); % set minimal value for generators rated power
Ert = mpc.bus(ssgen,VM).*baseV(ssgen); % [Vrms] machine no-load output voltage at frequency w=ws
J_GAIN_CNST = 0.005; % [s^2] constant defining moment of inertia
KD_GAIN_CNST = 0; 100; % constant defining damping factor
DROOP_R_GAIN_CNST = 20; % [rad/sec] constant defining droop control sloop parameter
% Other default parameters are based on this data:
Pm_rt = 3*Prt; % [W] rated mechanical (input) power
poles = 2; % number of machine poles (must be even)
J = J_GAIN_CNST*Prt/ws; % rotor moment of inertia
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
droop_Pref = mpc.gen(2:end,PG)*1e6; % [W] reference power for droop control (nominal generator power, single phase)
droop_R = DROOP_R_GAIN_CNST./Prt; % [rad/sec / W] droop control sloop parameter

% Create a new database for constructing the transmission network dynamic model
mpcC = line2shunt(mpc); % converts line charging capacitors to shunt elements
% Calculate capacitors to be placed at the sync machine outputs in simulation:
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

% Prepare network database for power flow calculations
mpcE = mpcC; % new database
% Add new buses to represent the synchronous machines behind
% their synchronous impedances.
% Define size and indices
nG = length(ssgen); % number of generators not including bus 1
[nRb, nCb]=size(mpcC.bus); % original size of bus matrix
[nRr, nCr]=size(mpcC.branch); % original size of branch matrix
[nRg, nCg]=size(mpcC.gen); % original size of generator matrix
indG = ((nRb+1):(nRb+nG))'; % indices of new buses
indR = ((nRr+1):(nRr+nG))'; % indices of new branches

% Copy original bus data
mpcE.bus = zeros(nRb + nG,nCb);
mpcE.bus(1:nRb,1:nCb) = mpcC.bus;
% Define new buses that represent the sync. machine voltage sources
mpcE.bus(indG,BUS_I)= indG; % new bus indices
mpcE.bus(indG,BUS_TYPE)= 2; % define new buses as PV buses
mpcE.bus(ssgen,BUS_TYPE)= 1; % define original buses as PQ buses
mpcE.bus(indG,PD)= 0; % real power demand
mpcE.bus(indG,QD)= 0; % reactive power demand
mpcE.bus(indG,BUS_AREA)= 1; 
mpcE.bus(indG,VM) = Ert./baseV(ssgen);
mpcE.bus(indG,BASE_KV)= baseV(ssgen)/1e3; 
mpcE.bus(indG,ZONE)= 1; 
mpcE.bus(indG,VMAX)= 1.2; 
mpcE.bus(indG,VMIN)= 0.8;
% Add new branches to connect the new buses above
% Copy original data
mpcE.branch = zeros(nRr+nG,nCr);
mpcE.branch(1:nRr,1:nCr) = mpcC.branch;
% Define new branches
mpcE.branch(indR,F_BUS) = indG; % new buses
mpcE.branch(indR,T_BUS) = ssgen; % original generator buses
mpcE.branch(indR,BR_R) = Ra./baseZ(ssgen);
mpcE.branch(indR,BR_X) = (Ls*ws)./baseZ(ssgen);
mpcE.branch(indR,BR_STATUS) = 1; %active line
mpcE.branch(indR,ANGMIN) = -360;
mpcE.branch(indR,ANGMAX) = +360;
% Define synchronous generators on the new buses
% Copy original data
mpcE.gen = zeros(nG+1,nCg);
mpcE.gen(1,:) = mpcC.gen(1,:); % copy infinite bus generator
mpcE.gen(2:(nG+1),GEN_BUS) = indG; % new generators on buses 15,16
mpcE.gen(2:(nG+1),PG) = droop_Pref  / 1e6;
mpcE.gen(:,QMAX) = 10 * max(droop_Pref  / 1e6);
mpcE.gen(:,QMIN) = -10 * max(droop_Pref  / 1e6);
mpcE.gen(2:(nG+1),VG) = mpcE.bus(indG,VM);
mpcE.gen(2:(nG+1),MBASE) = baseVA / 1e6;
mpcE.gen(2:(nG+1),GEN_STATUS) = 1;
mpcE.gen(2:(nG+1),PMAX) = Prt/1e6;
mpcE.gen(:,PMIN) = 0;
% Add shunt capacitors connected to syncronous generators
mpcE.bus(ssgen,BS) = mpcE.bus(ssgen,BS) +Cbus_ext.*baseV(ssgen).^2 *(ws/1e6);
% Represent loads as shunt elements instead of power sinks
mpcE.bus(:,PD) = 0;
mpcE.bus(:,QD) = 0;
Rtmpvec = [R_bus ; zeros(nG,1)];
Ltmpvec =  [L_bus ; inf*ones(nG,1)];
baseVtmpvec = [baseV ; baseV(ssgen) ];
mpcE.bus(:,GS) = mpcE.bus(:,GS) + (baseVtmpvec.^2).*real(1./(Rtmpvec-1j*ws*Ltmpvec))/1e6;
mpcE.bus(:,BS) = mpcE.bus(:,BS) - (baseVtmpvec.^2).*imag(1./(Rtmpvec-1j*ws*Ltmpvec))/1e6;
% Run power flow
mpcE = runpf(mpcE);
if (mpcE.success == 0)
    fprintf('%s\n\n',filename); % display network name    
    disp('power flow did not converge'); return;
end

% Compute initial capacitor voltages for simulation - using the power flow solution
vppa=(2^0.5)*mpcE.bus(ssgen,VM).*baseV(ssgen).*exp(1j*mpcE.bus(ssgen,VA)*pi/180);
vcap_init_d = real(vppa); % [V]
vcap_init_q = imag(vppa); % [V]
vcap_init_0 = 0*vppa; % [V]

% Signal routing for simulation
units_to_net = [1,4:(3+nG),2,(4+nG):(3+2*nG),3,(4+2*nG):(3+3*nG)];
[~,net_to_units] = sort(units_to_net);

% A short status report
fprintf('%s\n',filename);
fprintf('Number of buses = %d\n\n',N);
fprintf('Please wait ...\n');
fprintf('Simulation time = %d sec\n\n',round(sim_time));
if (N>100)
    fprintf('Simulation converging to steady state...\n');
    fprintf('Warning: long simulation time. This may take several minutes...\n');
elseif (N>20)
    fprintf('Simulation converging to steady state...\n');
end

% Run simulation (Simulink)
fprintf('Starting Simulink...\n\n');
exDroopA_sim;
set_param(bdroot,'StopTime',num2str(sim_time));
set_param(bdroot,'RelTol',num2str(sim_rel_tol));
set_param(bdroot,'MaxStep',num2str(sim_step));
% sim(bdroot);

% Small-signal analysis
% Simulate model and linearize at operating point:
io = getlinio(bdroot);
linsys = linearize(bdroot,io,sim_time);
[Ap,Bp,Cp,Dp]=ssdata(linsys); % linear system
% The linear model inputs and outputs are marked in Simulink.
% Input pertubation at droop-control referenece power
% Output measuement at transmission network dq0 currents
eeg = eig(Ap); % compute eigenvalues
fprintf('Eigenvalue maximal real part is %f\n',max(real(eeg)));


%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Visualize results:
%%%%%%%%%%%%%%%%%%%%%%%%
set(0,'defaulttextinterpreter','latex')
set(0,'defaultfigurecolor',[1 1 1])
set(0,'defaultaxesfontsize',9);

figure(1);
spy(Ap); 
title('state matrix A of linear system');
% Pole map
figure(2);
plot(real(eeg),imag(eeg),'kx'); xlim([-100 0]);
xlabel('Real Axis'); ylabel('Imaginary Axis');
title('Pole map');
% Bode plots
sys_inputs = 1;
sys_outputs = 1:2;
figure(3);
bode(linsys(sys_outputs,sys_inputs));


% Time-domain results
figure(1);
numplots = 5;
cur_plot = 1;
xmin = 0; % start time for plot
xmax = 6.5; sim_time; % final time for plot

plot(tsim,Pm/1e6,'-k','Color',[0 0 0],'LineWidth',0.7);
xlim([xmin xmax]);
ylabel('$P_m$ [MW]','FontSize',9);
xlabel('Time [s]');

axesHandles = findall(0,'type','axes');
set(axesHandles,'TickLabelInterpreter', 'latex')
