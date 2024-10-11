% Simulates the dynamics of a network that contains
% a long power line connecting to a photovoltaic source.
% The long transmission line and the photovoltaic source are modeled
% based on dq0 quantities.
%
% ** To start the simulation run this script.
% ** Associated Simulink file: 'exLongLinePV_sim'
%
% The network in this example includes two buses, which are
% connected by a long transmission line. The purpose of the simulation
% is to demonstrate the use of dq0 signals to model the photovoltaic
% source and its control, and to show the effects of the long power
% line (especially delays) on the stability of the system.
%
% bus 1 and bus 2 are connected by a long transmission line.
% bus 1 - infinite bus, modeled by a voltage source.
% bus 2 - renewable generator, modeled using dq0 signals.
% the model is based on paper "Y. Levron, S. Canaday, and
% R. W. Erickson, “Bus voltage control with zero
% distortion and high bandwidth for single-phase
% solar inverters,” IEEE Transactions on Power Electronics,
% 31 (1), pp. 258–269, Jan. 2016"

% The complete system is simulated in the time domain,
% using Simulink. The user may trigger a dynamic
% response by defining a step in input power.
%
% For more information about dq0 signals and modeling,
% refer to documentation in function 'ssNetw'

clc;

% Simulation parameters
sim_time = 1.5;  % [s] simulation time
sim_rel_tol = 1e-4;  % simulation accuracy
sim_step = sim_time/1e3; % simulation max step size
ws = 2*pi*50; % [rad/s] nominal system frequency
% Define step in input power:
% (to cancel the step define: rel_step=0)
rel_step = 0.02; % step size relative to the generator initial power
step_time = 1; % [s] when the step occurs.

% Infinite bus voltage:
bus1_vd = (2^0.5)*100e3; % [V] bus 1 (infinite bus) voltage, d component
bus1_vq = 0;  % [V] bus 1 (infinite bus) voltage, q component

% Long power line model
Rx = 4e-4; % [Ohm/m]
Lx = 1e-6; % [H/m]
Cx = 5e-12; % [F/m]
len = 200e3;  % [m]
Nn = 10; % number of sections

expected_delay_of_long_line = len*(Lx*Cx)^0.5;
% Transmission line state-space model
[At, Bt, Ct, Dt] = longLine(Rx, Lx, Cx, len, Nn, ws);
At=full(At); Bt=full(Bt);
Ct=full(Ct); Dt=full(Dt);

% Photovoltaic source model
Pren = 5000e3; % [W] average output power of photovoltaic source for a single phase (total power is 3 times higher)
Cgen = 200e-9; % [F] inverter output capacitance
Cbus = 1000e-6; % [F] bus capacitor (at the inverter input)
Vbus_ref = 800; % [V] bus voltage reference
kp = 1e-5;  % proportional constant of inverter PI controller 
ki = 1e-5; % integral constant of inverter PI controller

% Run simulation (Simulink)
exLongLinePV_sim;
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

figure;
numplots = 3;
cur_plot = 1;
xmin = 0.8;
xmax = 1.5;

subplot(numplots,1,cur_plot);
plot(tsim,Pdc/1e3,'-k','Color',[0 0 0],'LineWidth',0.5);
xlim([xmin xmax]);
ylabel('$P_{pv}$ [kW]','FontSize',9)
cur_plot = cur_plot+1;

subplot(numplots,1,cur_plot);
plot(tsim,genp(:,1),'-k','Color',[0 0 0],'LineWidth',0.5);
xlim([xmin xmax]);
ylabel('$P_2$ [kW]','FontSize',9)
cur_plot = cur_plot+1;

subplot(numplots,1,cur_plot);
plot(tsim,gen_cur_dq0(:,1),'-k','Color',[0 0 0],'LineWidth',0.5);
xlim([xmin xmax]);
xlabel('Time [s]','FontSize',9)
ylabel('$I_{d,2}$ [A]','FontSize',9)

axesHandles = findall(0,'type','axes');
set(axesHandles,'TickLabelInterpreter', 'latex')
