% Stability analysis of a 2-bus power system, using dq0 signals.
% This file demonstrates:
% -- How to use the function 'closedLoop' to construct a
% state-space model of a complete power system
% (transmission network + generators/loads).
% -- How to test stability theoretically by examining eigenvalues of
% the state matrix A.
%
% ** To start the simulation - run this script.
% ** Associated simulink file: 'ex2busStability_sim'
%
% The user may choose a scenario to simulate:
% - Scenario 1: stable system
%       In this case the system converges to its operating point.
%       System dynamics may be examined by stepping the
%       mechnical input power.
% - Scenario 2: unstable system
%       In this case the system is unstable, and enters a limit cycle
%       around its operating point.
%
% Stability is controlled by tuning the branch resistance R.
% Higher resistance means higher energy losses, and better stability.
%
% The network in this example includes two buses:
% bus 1 - infinite bus, modeled by a voltage source.
% bus 2 - synchronous machine (generator)
% These are connected through a single branch
% with inductance L and resistance R.

clc;
close all;

% Choose_scenario
% 1 = stable system
% 2 = unstable system
scenario = 2;

switch scenario
    case 1 % stable system, with input step
        R = 1; % [Ohm] network branch resistance
            % (with this resistance the system is stable)
        L = 0.5; % [H] network branch inductance
        step_time = 9; % [s] time of step in mechanical power
        step_size = 30e3; % [W] step size
    case 2 % unstable system (entering a limit cycle)
        R = 0; % [Ohm] network branch resistance
            % (the resistance is too low and the system is unstable)
        L = 0.5; % [H] network branch inductance
        step_time = 9; % no input step
        step_size = 30e3; 
    otherwise
        display('Unknown scenario');
        return;
end


% System parameters
% Operating point parameters:
V1mag = 10e3; % [Vrms] voltage magnitude (RMS) of bus 1 (infinite bus)
Vgmag = 1.02 * V1mag; % [Vrms] voltage magnitude (RMS) of synchronous machine in bus 2
Pg = 200e3; % [W] active power (each phase) of synchronous machine in bus 2

% Simulation parameters:
sim_time = 10;  % [s] simulation time
sim_rel_tol = 1e-6;  % simulation accuracy
sim_step = sim_time/1e4; % simulation max step size
ws = 2*pi*50; % [rad/s] nominal system frequency

% Compute the system operating point
% The power flow is computed analytically
% (in more complex networks the power flow equations are usually solved numerically)
% Define constants and solve quadratic equation:
a1 = Vgmag^2;
a2 = R*V1mag;
a3 = ws*L*V1mag;
a4 = ((ws*L)^2+R^2)*Pg - R*Vgmag^2;
vopy = (a3*a4 + a2*(a1*a2^2 + a1*a3^2 - a4^2)^(1/2))/(a2^2 + a3^2);
vopx = (Vgmag^2 - vopy^2)^0.5;
iopx = ((vopx-V1mag)*R + ws*L*vopy)/((ws*L)^2+R^2);
iopy = (vopy*R+ws*L*(V1mag-vopx))/((ws*L)^2+R^2);
Qg_op = vopy*iopx - vopx*iopy; % reactive power of synchronous machine in bus 2
ph2_op = atan2(vopy,vopx); % phase of synchronous machine in bus 2
ph2_deg = ph2_op*180/pi; % same phase in degrees
% Operating point for Simulink:
vol_dq0_op = 2^0.5*[V1mag; vopx; 0; vopy; 0; 0]; % operating point of voltage dq0 signals
cur_dq0_op = 2^0.5*[-iopx; iopx; -iopy; iopy; 0; 0]; % operating point of current dq0 signals

% Synchronous machine parameters
Pg_max = 12*Pg;  % [W] maximal generator active power output (total of three-phases)
Pm_op = 3*Pg; % [W] mechanical input power at operating point
poles = 2; % number of machine poles
lamda = (2^0.5) * Vgmag/ws; % machine voltage gain (voltage / angular velocity)
J_GAIN_CNST = 0.001;
J = J_GAIN_CNST*Pg_max/ws; % rotor moment of inertia
KD_GAIN_CNST = 10;
Kd = KD_GAIN_CNST*J*ws; % damping factor
poles_2Jws = poles./(2*ws*J); % constant used in simulink

% Construct linear (small-signal) model of infinite bus (bus 1)
[A1, B1, C1, D1, G1] = infBus();

% Construct linear (small-signal) model of synchronous machine (bus 2)
params = struct([]);
params(1).poles = poles;
params.J = J;
params.Kd = Kd;
op = struct([]);
op(1).Vmag = Vgmag;
op.ph = ph2_op;
op.P = Pg;
op.Q = Qg_op;
[A2, B2, C2, D2, G2] = syncMachA(params, op, ws);

% Construct linear model of the transmission network
% (the network in this example is a single branch
% with inductance L and resistance R)
R_bus = [0; 0];
L_bus = [0; 0];
C_bus = [0; 0];
Rtil_bus = [0; 0];
Tau = zeros(2);
Lb = [0 L; 0 0];
Rb = [0 R; 0 0];
[Anet, Bnet, Cnet, Dnet] = ssNetw(R_bus, L_bus, C_bus, Rtil_bus, Rb, Lb, Tau, ws);

% Construct linear (small-signal) model of the
% complete network (generator + network)
unitData = struct([]);
unitData{1}.A=A1; unitData{1}.B=B1; unitData{1}.G=G1;
unitData{1}.C=C1; unitData{1}.D=D1;
unitData{2}.A=A2; unitData{2}.B=B2; unitData{2}.G=G2;
unitData{2}.C=C2; unitData{2}.D=D2;
% Construct state-space model:
[A, B, C, D] = closedLoop(Anet, Bnet, Cnet, Dnet, unitData);

% Test stability of the power system by
% examining the eigenvalues of the state matrix A
dr_max = eigs(A,1,0.1);
% Note: 'eigs' is prefered over 'eig' when A is large and sparse
fprintf('Largest eigenvalue real part is %f\n',real(dr_max));
if (real(dr_max) < 0)
    disp('System is stable at operating point.');
else
    disp('System is unstable at operating point!');
end

% Run simulation (Simulink)
ex2busStability_sim;
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
tmax = sim_time;
figure(1);
switch scenario
    case 1, Nplot = 4; tmin = 8;
    case 2, Nplot = 2; tmin = 9.8;
end

subplot(Nplot,1,1);
plot(tsim,genp*1e3,'-k','Color',[0 0 0],'LineWidth',0.5);
xlim([tmin tmax]);
ylabel('$P_1,P_2$ [kW]');
% Compare results to powers computed by linear model:
genplin = 0*genp;
genplin(:,1) = 0.5*(lin_v(:,1).*lin_i(:,1)+lin_v(:,3).*lin_i(:,3)); % active power bus 1
genplin(:,2) = 0.5*(lin_v(:,2).*lin_i(:,2)+lin_v(:,4).*lin_i(:,4)); % active power bus 2
genplin(:,1) = genplin(:,1); % convert to kW
genplin(:,2) = genplin(:,2);
hold on;
plot(tsim,genplin,'--k','Color',[0 0 0],'LineWidth',0.8);
hold off;

switch scenario
    case 1, title('Stable system - response to step in mechanical power');
    case 2, title('System is unstable and enters a limit cycle around its operating point');
end

subplot(Nplot,1,2);
plot(tsim,genq,'-k','Color',[0 0 0],'LineWidth',0.5);
xlim([tmin tmax]);
ylabel('$Q_1,Q_2$ [kVAr]');
% Compare results to powers computed by linear model:
genqlin = 0*genp;
genqlin(:,1) = 0.5*(lin_v(:,3).*lin_i(:,1)-lin_v(:,1).*lin_i(:,3)); % active power bus 1
genqlin(:,2) = 0.5*(lin_v(:,4).*lin_i(:,2)-lin_v(:,2).*lin_i(:,4)); % active power bus 2
genqlin(:,1) = genqlin(:,1) / 1e3; % convert to kVAr
genqlin(:,2) = genqlin(:,2) / 1e3;
hold on;
plot(tsim,genqlin,'--k','Color',[0 0 0],'LineWidth',0.86);
hold off;

if (scenario==1)
% Current signals in the abc reference frame:
    subplot(Nplot,1,3);
    plot(tsim,gen_cur_abc,'LineWidth',0.5);
    xlim([tmin tmax]);
    ylabel('$I_{abc,2}$ [A]');
    
% Normalized error in linear model voltage signals:
    subplot(Nplot,1,4);
    plot(tsim,voltage_error./repmat(abs(vol_dq0_op'),length(tsim),1),'LineWidth',0.5);
    xlim([tmin tmax]);
    ylabel('Error');
end

xlabel('Time [s]');
axesHandles = findall(0,'type','axes');
set(axesHandles,'TickLabelInterpreter', 'latex')
