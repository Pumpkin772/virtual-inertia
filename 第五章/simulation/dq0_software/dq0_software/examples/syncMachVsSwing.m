% Comparing the transient response of the physical
% and approximate (swing equation) synchronous machine models.
% The transient is triggered with a step in mechanical (input) power.
% Both models are tested when connected to an infinite bus.
%
% ** Associated Simulink file: 'syncMachVsSwing_sim'
%
% The physical (accurate) model is based on the
% book "Electric Machinery" by Fitzgerald. The presented
% model is based on dq0 quantities, and is identical to the
% model presented in the book, except that the input
% voltages and output currents are represented in
% a references frame rotating with angle ws*t,
% (instead of theta - the rotor electrical angle),
% where ws is the frequency of the infinite bus.
% This change of coordinates allows direct connection of
% the model to the infinite bus, or to any other model
% (network / load / generator) that is represented in
% the same reference frame.
%
% The approximated model is based on the (classical)
% swing equation, and includes the effects of a series
% synchronous inductance and armature resistance.
% Similar to the physical model, the input
% voltages and output currents are represented in
% a references frame rotating with angle ws*t.
%
% Main results:
% the two models are identical at steady-state,
% but have different dynamics, and respond
% differently during transients.

% Simulation parameters
sim_time = 10;  % [s] simulation time
sim_rel_tol = 1e-4;  % simulation accuracy
sim_step = sim_time/1e4; % simulation max step size

% Define step in input power:
% (to cancel the step define: rel_step=0)
step_time = 5; % [s] when the step occurs
rel_step = 1; % relative step size in mechanical power

% Power system parameters
ws = 2*pi*50; % [rad/s] nominal system frequency (infinite bus frequency)

% Synchronous machine(s) parameters
% Enter the following specifications:
Prt = 10e6; % [W] machine rated power (single phase) (maximum power for normal operation)
Ert = 10e3; % [Vrms] machine rated voltage (RMS)
% Default parameters based on this data:
Pm_rt = 3*Prt; % [W] rated mechanical (input) power
poles = 2; % number of machine poles (must be even)
J_GAIN_CNST = 0.005; % [s^2] constant defining moment of inertia
KD_GAIN_CNST = 10; % constant defining damping factor
J = J_GAIN_CNST*Prt/ws; % rotor moment of inertia
Kd = KD_GAIN_CNST*J*ws; % damping factor
poles_2Jws = poles/(2*J*ws); % constant used in simulation
If_dc = 0.2*Prt/Ert; % [A] field winding nominal DC current
Lpu = Ert^2/(Prt*ws); % [H] per-unit inductance
Ld = 0.6*Lpu; % [H] direct axis synchronous inductance
Lq = Ld; % [H] quadrature axis synchronous inductance
L0 = 0.1*Ld; % [H] zero sequence inductance
Ra = 0.05*Ert^2/Prt; % [Ohm] armature resistance
Laf = (2^0.5)*Ert/(ws*If_dc); % [H] stator to rotor mutual inductance
Lff = 2*Laf^2/Ld; % [H] field winding self inductance
Rf = Lff/1.2; % [Ohm] field winding resistance
Vf_dc = If_dc*Rf; % [V] field winding nominal DC voltage
cc1 = 2*Ld*Lff - 3*Laf^2; % constant used in simulation
Ls = (Ld+Lq)/2; % [H] synchronous inductance (used in simplified sync. machine model)

% Run simulation (Simulink)
syncMachVsSwing_sim;
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
numplots = 5;
cur_plot = 1;
xmin = 4.5; % start time for plot
xmax = 6.5; % final time for plot

subplot(numplots,1,cur_plot);
plot(tsim,p1/1e6,'-k','Color',[0 0 0],'LineWidth',0.7);
hold on;
plot(tsim,p2/1e6,'--k','Color',[0.15 0.15 0.15],'LineWidth',0.7);
xlim([xmin xmax]);
ylabel('$P$ [MW]');
title('solid - physical model; dashed - approximate model');
cur_plot = cur_plot+1;

subplot(numplots,1,cur_plot);
plot(tsim,f1,'-k','Color',[0 0 0],'LineWidth',0.7);
hold on;
plot(tsim,f2,'--k','Color',[0.15 0.15 0.15],'LineWidth',0.7);
xlim([xmin xmax]);
ylabel('$f$ [Hz]');
cur_plot = cur_plot+1;

subplot(numplots,1,cur_plot);
plot(tsim,delta1,'-k','Color',[0 0 0],'LineWidth',0.7);
hold on;
plot(tsim,delta2,'--k','Color',[0.15 0.15 0.15],'LineWidth',0.7);
xlim([xmin xmax]);
ylabel('$\delta$ [deg]');
cur_plot = cur_plot+1;

subplot(numplots,1,cur_plot);
plot(tsim,idq0_1(:,[1 2]),'-k','Color',[0 0 0],'LineWidth',0.7);
hold on;
plot(tsim,idq0_2(:,[1 2]),'--k','Color',[0.15 0.15 0.15],'LineWidth',0.7);
xlim([xmin xmax]);
ylabel('$I_{dq}$ [A]');
ylim([-400 2200]);
cur_plot = cur_plot+1;

subplot(numplots,1,cur_plot);
plot(tsim,pm/1e6,'-k','Color',[0 0 0],'LineWidth',0.7);
xlim([xmin xmax]);
ylabel('$P_m$ [MW]');
ylim([10 35]);
cur_plot = cur_plot+1;

xlabel('Time [s]');

axesHandles = findall(0,'type','axes');
set(axesHandles,'TickLabelInterpreter', 'latex')
