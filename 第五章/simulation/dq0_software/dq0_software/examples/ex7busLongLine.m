% Simulates the dynamics of a 7-bus system with long transmission lines,
% synchronous generators, and photovoltaic generators. All models are
% based on dq0 quantities.
% The simulation demonstrates the use of dq0 signals to model
% all the system components, and shows the effects of the long power
% line (especially delays) on the dynamics and stability of the system.
%
% This example also demonstrates how to model a photovoltaic generator
% based on virtual inertia. This generator uses a local energy storage to
% emulate the inertia of a spinning machine, with an objective to
% stabilize the system. It is shown that if the virtual moment of inertia
% is too low, the system is unstable. To see this result, set a low value
% to the parameter JPV_CNST (for instance JPV_CNST =0.0001)
%
% ** Associated Simulink file: 'ex7busLongLine_sim'
%
% SYSTEM STRUCTURE
% Buses
% bus 1 - infinite bus (constant frequency and voltage)
% bus 2 - load
% bus 3 - load
% bus 4 - synchronous generator, physical model, droop control, connected through a transformer
% bus 5 - photovoltaic generator, control based on "virtual inertia" (see below),  connected through a transformer
% bus 6 - load
% bus 7 - load
%
% "virtual inertia" control:  See paper "Synchronverters: Inverters That Mimic
% Synchronous Generators" and references within.
%
% Lines (all lines are modeled as long transmission lines)
% 1 -> 2
% 1 -> 4
% 2 -> 3
% 3 -> 4
% 2 -> 5
% 4 -> 6
% 6 -> 7


clc;

% Simulation parameters
sim_time = 5;  % [s] simulation time
sim_rel_tol = 1e-7;  % simulation accuracy
sim_step = sim_time/1e4; % simulation max step size

% Transient simulation: define input steps
tss2 = 4.5; % switch time for extra load on bus 2 (to cancel define as 1e9)
tss3 = 1e9; % switch time for extra load on bus 3 (to cancel define as 1e9)
tss6 = 1e9; % switch time for extra load on bus 6 (to cancel define as 1e9)
tss7 = 1e9; % switch time for extra load on bus 7 (to cancel define as 1e9)
tss4 = 1e9; % step time for power set-point of synchronous generator on bus 4 (to cancel define as 1e9)
rel_step_4 = 0.5; % relative step size for power set-point of synchronous generator
tss5 = 4; % step time for power of photovoltaic generator on bus 5 (to cancel define as 1e9)
rel_step_5 = -0.5; % relative step size for power of photovoltaic generator

% System and transmission lines parameters
ws = 2*pi*50; % [rad/s] nominal system frequency
Pbase = 300e6; % [W] per-unit base, and the rated power of the transmission lines.
Vbase = 161e3; % [Vrms] per-unit base, and the rated voltage of the transmission lines.
Lx = 1.35e-6; % [H/m] transmission lines inductance per unit length
Cx = 8.45e-12; % [F/m] transmission lines capacitance per unit length
Rx = 1.07e-4; % [Ohm/m] transmission lines resistance per unit length

% Transmission line parameters
% lenXY - length of line from bus X to bus Y
len12 = 227e3; % [m] line length
NN12 = 5; % number of sections in line model
len14 = 302e3; % [m] line length
NN14 = 5; % number of sections in line model
len23 = 54e3; % [m] line length
NN23 = 3; % number of sections in line model
len34 = 47e3; % [m] line length
NN34 = 2; % number of sections in line model
len25 = 85e3; % [m] line length
NN25 = 3; % number of sections in line model
len46 = 78e3; % [m] line length
NN46 = 3; % number of sections in line model
len67 = 63e3; % [m] line length
NN67 = 2; % number of sections in line model

% Unit (bus) parameters
% bus 1 - infinite bus
vinf_d = 1.1*(2^0.5)*Vbase; % [V] voltage source d-component
vinf_q = 0; % [V] voltage source q-component
vinf_0 = 0; % [V] voltage source 0-component
% bus 2 - load
PL2 = -0.34*Pbase; % load active power (single phase)
QL2 = -0.24*Pbase; % load reactive power
C2 = 0.1*(len12/NN12+len23/NN23+len25/NN25)*Cx; % shunt capacitance on bus
% bus 3 - load
PL3 = -0.28*Pbase; % load active power (single phase)
QL3 = -0.08*Pbase; % load reactive power
C3 = 0.1*(len23/NN23+len34/NN34)*Cx; % shunt capacitance on bus
% bus 4 - synchronous generator, physical model (Fitzgerald), droop control, connected through a transformer
Prt = 0.5*Pbase; % [W] machine rated power (single phase) (maximum power for normal operation)
Ert = 1.1*22e3; % [Vrms] machine no-load output voltage at frequency w=ws
J_GAIN_CNST = 0.005; % [s^2] constant defining moment of inertia
KD_GAIN_CNST = 0; 100; % constant defining damping factor
DROOP_R_GAIN_CNST = 20; % [rad/s] constant defining droop control sloop parameter
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
% Synchronous generator droop control paramters
droop_Pref =  0.25*Pbase; % [W] reference power for droop control (nominal generator power, single phase)
droop_R = DROOP_R_GAIN_CNST./Prt; % [rad/sec / W] droop control sloop parameter
% Additional parameters for bus 4
tratio4 = 22e3 / Vbase; % transformer ratio
C4 = 0.1*(len14/NN14+len34/NN34+len46/NN46)*Cx; % shunt capacitance on bus
% bus 5 - photovoltaic generator, uses virtual inertia, connected through a transformer
Ppvrt = 0.2*Pbase; % photovoltaic generator rated power (single phase) (maximum power for normal operation)
PG5 = 0.8*Ppvrt; % nominal photovoltaic generator power, (single phase)
Epv = 22e3; % [Vrms]  photovoltaic generator no-load output voltage at frequency w=ws
JPV_CNST = 0.0005; % [s^2] constant defining the virtual moment of inertia
KDPV_CNST = 200; % [s^2] constant defining the damping factor
Jpv = JPV_CNST*Ppvrt/ws; % [kg*m^2] virtual moment of inertia (emulated by control)
Kdpv = KDPV_CNST*Jpv*ws;
Lspv = 0.1*Epv^2/(Ppvrt*ws); % [H] photovoltaic generator series inductance
Rspv =  0.005*Epv.^2./Ppvrt; % [Ohm] photovoltaic generator series resistance
% Additional parameters for bus 5
tratio5 = 22e3 / Vbase; % transformer ratio
C5 = 0.1*(len25/NN25)*Cx; % shunt capacitance on bus
% bus 6 - load
PL6 = -0.21*Pbase; % load active power (single phase)
QL6 = -0.04*Pbase; % load reactive power
C6 = 0.1*(len46/NN46+len67/NN67)*Cx; % shunt capacitance on bus
% bus 7 - load
PL7 = -0.07*Pbase; % load active power (single phase)
QL7 = -0.03*Pbase; % load reactive power
C7 = 0.1*(len67/NN67)*Cx; % shunt capacitance on bus

% Build transmission lines state-space models
[At12, Bt12, Ct12, Dt12] = longLine(Rx, Lx, Cx, len12, NN12, ws);
[At14, Bt14, Ct14, Dt14] = longLine(Rx, Lx, Cx, len14, NN14, ws);
[At23, Bt23, Ct23, Dt23] = longLine(Rx, Lx, Cx, len23, NN23, ws);
[At34, Bt34, Ct34, Dt34] = longLine(Rx, Lx, Cx, len34, NN34, ws);
[At25, Bt25, Ct25, Dt25] = longLine(Rx, Lx, Cx, len25, NN25, ws);
[At46, Bt46, Ct46, Dt46] = longLine(Rx, Lx, Cx, len46, NN46, ws);
[At67, Bt67, Ct67, Dt67] = longLine(Rx, Lx, Cx, len67, NN67, ws);

% Model loads as series R-L impedances
R2 = -PL2*Vbase^2/(PL2^2+QL2^2); L2 = -(1/ws)*QL2*Vbase^2/(PL2^2+QL2^2);
R3 = -PL3*Vbase^2/(PL3^2+QL3^2); L3 = -(1/ws)*QL3*Vbase^2/(PL3^2+QL3^2);
R6 = -PL6*Vbase^2/(PL6^2+QL6^2); L6 = -(1/ws)*QL6*Vbase^2/(PL6^2+QL6^2);
R7 = -PL7*Vbase^2/(PL7^2+QL7^2); L7 = -(1/ws)*QL7*Vbase^2/(PL7^2+QL7^2);
[AL2, BL2, CL2, DL2] = balLoadsRL(R2, L2, ws); % series RL load
[AL3, BL3, CL3, DL3] = balLoadsRL(R3, L3, ws); % series RL load
[AL6, BL6, CL6, DL6] = balLoadsRL(R6, L6, ws); % series RL load
[AL7, BL7, CL7, DL7] = balLoadsRL(R7, L7, ws); % series RL load
CL2=-CL2; DL2=-DL2; % reverse the direction of current - positive when injected into network
CL3=-CL3; DL3=-DL3; % reverse the direction of current - positive when injected into network
CL6=-CL6; DL6=-DL6; % reverse the direction of current - positive when injected into network
CL7=-CL7; DL7=-DL7; % reverse the direction of current - positive when injected into network

%%%%%%%%%%%% run simulation (Simulink) %%%%%%%%%%
fprintf('Starting Simulink... please wait... \n\n');
ex7busLongLine_sim;
set_param(bdroot,'StopTime',num2str(sim_time));
set_param(bdroot,'RelTol',num2str(sim_rel_tol));
set_param(bdroot,'MaxStep',num2str(sim_step));
tic
sim(bdroot);
toc


%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Visualize results:
%%%%%%%%%%%%%%%%%%%%%%%%
set(0,'defaulttextinterpreter','latex')
set(0,'defaultfigurecolor',[1 1 1])
set(0,'defaultaxesfontsize',9);

figure(1);
numplots = 5;
cur_plot = 1;
xmin = 3.8; % step_time-0.5; % start time for plot
xmax = tsim(end); % step_time+1; % final time for plot

subplot(numplots,1,cur_plot);
plot(tsim,[f4 f5],'Color',[0 0 0],'LineWidth',0.6);
xlim([xmin xmax]);
ylabel('$f$ [Hz]','FontSize',9);
cur_plot = cur_plot+1;

subplot(numplots,1,cur_plot);
plot(tsim,p3/1e6,'Color',[0 0 0],'LineWidth',0.6);
xlim([xmin xmax]);
ylabel('$P_3$ [MW]','FontSize',9);
cur_plot = cur_plot+1;

subplot(numplots,1,cur_plot);
plot(tsim,p4/1e6,'Color',[0 0 0],'LineWidth',0.6);
xlim([xmin xmax]);
ylabel('$P_4$ [MW]','FontSize',9);
cur_plot = cur_plot+1;

subplot(numplots,1,cur_plot);
plot(tsim,p5/1e6,'Color',[0 0 0],'LineWidth',0.6);
xlim([xmin xmax]);
ylabel('$P_5$ [MW]','FontSize',9);
cur_plot = cur_plot+1;

xlabel('Time [s]');

axesHandles = findall(0,'type','axes');
set(axesHandles,'TickLabelInterpreter', 'latex')
