% Modeling and simulation of a three-phase
% photovoltaic inverter based on dq0 signals.

% ** Associated Simulink file: 'pv_dq_sim'

% Note: the inverter model is based on averaged signals.
% - All signals are averaged in respect to the high switching
% frequency of the power devices.
% - All signals are represented in the dq0 domain.

% Modeling assumptions:
% ---------------------
% 1. Typical dq0-based PQ control scheme similar to the one
% described in N. Kroutikova, C. A. Hernandez-Aramburo, and T.
% C. Green, “Statespace model of grid-connected inverters under
% current control mode,” IET Electr. Power App., vol. 1, no. 3,
% pp. 329–338, 2007.
% 2. Standard average-signal model. High freuqency 
% switching harmonics are ignored.
% 3. The inverter is connected to an infinite bus,
% with a constant frequency ws.
% 4. All three-phase voltages and currents are balanced (v0=0, i0=0).
% 5. The PLL loop is ideal and perfectly extracts the 
% dq0 components and frequency of the grid voltage with
% zero distortion.
% 6. the MPPT unit and DC-DC converter are ideal an
% represented by a single constant power source (pdc).
% 7. All power conversion stages are assumed ideal
% (lossless and no internal energy storage).
% 8. Bus capacitor dynamics are approximated
% by a linear model: d(Vdc)/dt = (power charging cap.)/(C*Vdcset).
% 9. inverter AC output voltage is approximated
% by v_d = Vdcset*d_d , v_q = Vdcset*d_q, 
% ignoring variations in the bus capacitor voltage.
% 10. Current controler uses common "cross" topology,
% as described in the paper mentioned above.
% 11. The current controller and active-power compensators
% are based on simple proportional-integral controllers.

clc;
clear variables;

% Simulation parameters:
sim_time = 2;  % [s] simulation time
sim_rel_tol = 1e-4;  % simulation accuracy
sim_step = sim_time/1e4; % simulation max step size
% input ramp in DC power:
%(to cancel the ramp select ramp_size=0)
ramp_duration = 1; % [sec]
ramp_size = -0.1; % relative to DC power at steady-state
ramp_starts_at = 0.2; % [sec]

% Display: plot range:
tminplot = 0; % start time for plot
tmaxplot = sim_time; % final time for plot

% Global system parameters:
ws = 2*pi*50; % [rad/sec] nominal grid frequency
Vg = 22e3 / 3^0.5; % [Vrms] nominal grid voltage (line-to-ground RMS voltage)

% Photovoltaic inverter parameters
% Rated values:
Pinv_rt = 2e6; % [W] inverter rated power (total for three phase)
Vinv_rt = Vg; % [V] inverter rated output voltage (line-to-ground RMS voltage)
% The following parameters are assigned default values that represent
% a possible working design. All parameters are calculated based on
% the rated power and voltage, and will scale properly for different rated values.
% These values can be used "as is" for preliminary system testings,
% or can be modified for a more accurate analysis.
Qstar = 0; % [VAr] desired reactive power (per-phase)
L = 10*(1/ws)*Pinv_rt/Vinv_rt^2; % [H] total inductance (per phase) at inverter output
Ebus = 0.4e-3*Pinv_rt; % [J] total bus capacitor stored energy at nominal bus voltage
duty_max = 0.9; % theoretical maximal duty-cycle required to generate the nominal line-to-ground voltage at each phase.
% current loop compensator constants:
Kcp = 10*4*pi*L; % proportional (P) constant
Kci = Kcp^2/(2*L); % integral (I) constant
% active power loop compensator constants:
Kpp = 200*(Ebus/Vinv_rt); % proportional (P) constant
Kpi = 0.1*Kpp^2/(Ebus/Vinv_rt); % integral (I) constant
% constants for inverter model:
Vdcset = (2^0.5)*Vinv_rt/duty_max; % [V] desired bus voltage (transformed to bridge input)
C = 2*Ebus/Vdcset^2; % [F] effective bus capacitor (transformed to bridge input)

% system operating point data
Pdc_op = 0.8*Pinv_rt; % [W] inverter DC input power at steady-state (equals to the total three-phase output power)

% Run simulation (Simulink)
disp('Starting Simulink...')
pv_dq_sim;
set_param(bdroot,'StopTime',num2str(sim_time));
set_param(bdroot,'RelTol',num2str(sim_rel_tol));
set_param(bdroot,'MaxStep',num2str(sim_step));
disp('Transient simulation... Please wait ...')
sim(bdroot);
disp('Simulation completed.')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% plot results %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

set(0,'defaulttextinterpreter','latex')
set(0,'defaultfigurecolor',[1 1 1])
set(0,'defaultaxesfontsize',9);
cur_plot = 1;

numplots = 3;
figure(1);

subplot(numplots,1,cur_plot);
plot(tsim,pdc/1e3,'--k','Color',[0 0 0],'LineWidth',0.7); hold on;
plot(tsim,pgrid/1e3,'-k','Color',[0 0 0],'LineWidth',0.7); hold off;
xlim([tminplot tmaxplot]);
ylabel('$P$ [kW]','FontSize',12);
cur_plot = cur_plot+1;

subplot(numplots,1,cur_plot);
plot(tsim,pdc/1e3,'--k','Color',[0 0 0],'LineWidth',0.7); hold on;
plot(tsim,pgrid/1e3,'-k','Color',[0 0 0],'LineWidth',0.7); hold off;
xlim([tminplot tmaxplot]);
ylim([pdc(end)/1e3*0.999 pdc(end)/1e3*1.001]);
ylabel('$P$ [kW]','FontSize',12);
cur_plot = cur_plot+1;

subplot(numplots,1,cur_plot);
plot(tsim,tsim.^0,'--k','Color',[0 0 0],'LineWidth',0.7); hold on;
plot(tsim,vdc/Vdcset,'-k','Color',[0 0 0],'LineWidth',0.7); hold off;
xlim([tminplot tmaxplot]);
ylabel('$v_{dc}/V_{dc,set}$','FontSize',12);
ylim([0.9 1.05]);
cur_plot = cur_plot+1;

xlabel('Time [s]');
axesHandles = findall(0,'type','axes');
set(axesHandles,'TickLabelInterpreter', 'latex')

% additional text messages
disp(' '); disp('Additional information:');
disp(strcat(['Cbus (rated at 900V) is ',num2str(round((2*Ebus/900^2/1e-3)*100)/100)],' mF'));

