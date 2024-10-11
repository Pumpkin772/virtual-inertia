% Simulates the dynamics of a network containing
% a renewable energy source, using dq0 signals.
%
% ** To start the simulation - run this script.
% ** Associated Simulink file: 'ex3busPV_sim'
%
% The network in this example includes three buses:
% bus 1 - infinite bus, modeled by a voltage source
% bus 2 - renewable generator. This generator is
% modeled by a power source in parallel to a capacitor
% bus 3 - balanced three-phase resistive load
%
% The network includes two branches (power lines)
% branch 1->2: includes an inductor and a resistor
% branch 2->3: includes an inductor and a resistor
%
% The complete system is simulated in the time domain,
% using Simulink. The user may trigger a dynamic
% response by defining a step in input power.
%
% For more information about dq0 signals and modeling,
% refer to documentation in function 'ssNetw'

clc;

% Simulation parameters
sim_time = 0.1;  % [s] simulation time
sim_rel_tol = 1e-4;  % simulation accuracy
sim_step = sim_time/1e3; % % simulation max step size
ws = 2*pi*50; % [rad/s] nominal system frequency
% Define step in input power:
% (to cancel the step define: rel_step=0)
rel_step = 0.6; % step size relative to the generator initial power
step_time = 0.05; % [s] when the step occurs
% Bus parameters:
bus1_vd = 2^0.5*1e3; % [V] bus 1 (infinite bus) voltage, d component
bus1_vq = 0;  % [V] bus 1 (infinite bus) voltage, q component
Pren = 0.5e3; % [W] bus 2 initial (per-phase) output power of renewable source
gen_vd_init = bus1_vd; % [V] bus 2 generator output capacitor initial voltage, d component
gen_vq_init = bus1_vq; % [V] bus 2 generator output capacitor initial voltage, q component
Cgen = 2e-3*Pren/(gen_vd_init^2+gen_vq_init^2); % [F] bus 2 output capacitor of renewable source
Rload = 1e3; % [Ohm] bus 3 load resistance

% Build dynamic model of the network
% (network in this example includes two branches: 1->2 and 2->3)
R12 = 10; % [ohm] branch 1->2 resistance
R23 = 1; % [ohm] branch 2->3 resistance
L12 = 5; % [H] branch 1->2 inductance
L23 = 0.06; % [H] branch 2->3 inductance
R_bus = [0 0 0];
L_bus = [0 0 0];
C_bus = [0 0 0];
Rtil_bus = [0 0 0];
Lb = [0 L12 0; 0 0 L23; 0 0 0]; 
Rb = [0 R12 0; 0 0 R23; 0 0 0];
Tau = sparse(3,3);
% Build state-space model of the network in dq0 coordinates:
[A, B, C, D] = ssNetw(R_bus, L_bus, C_bus, Rtil_bus, Rb, Lb, Tau, ws);
A=full(A); B=full(B);
C=full(C); D=full(D);

% Run simulation (Simulink)
ex3busPV_sim;
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
Nplot = 3;

subplot(Nplot,1,1);
plot(tsim,genp,'LineWidth',0.5);
ylabel('$P(t)$ [kW]');

subplot(Nplot,1,2);
plot(tsim,genq,'LineWidth',0.5);
ylabel('$Q(t)$ [kVAr]');

subplot(Nplot,1,3);
plot(tsim,gen_cur_abc,'LineWidth',0.5);
ylabel('$I_{abc,2}$ [A]');

xlabel('Time [s]');

axesHandles = findall(0,'type','axes');
set(axesHandles,'TickLabelInterpreter', 'latex')
