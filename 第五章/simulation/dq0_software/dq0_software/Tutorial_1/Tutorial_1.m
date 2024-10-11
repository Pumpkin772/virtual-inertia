% This script is used as part of Tutorial 1.
% Detailed explanations may be found in the manual.

close all;

% Load network data
load sdq0_inputs.mat; % this file is created by the graphical-user-intertface

ws = 2*pi*50; % [rad/s] nominal system frequency 
N = size(Lb,1);  % number of buses in the network

%% Define parameters for Simulink model
% Create dq0 state-space model of the transmission network
[A1,B1,C1,D1] = ssNetw(R_bus,L_bus,C_bus,Rtil_bus,Rb,Lb,Tau,ws); % build state-space model without R-L loads
[A2,B2,C2,D2] = balLoadsRL(R_bus_SLD,L_bus_SLD,ws); % build state-space model of R-L loads
[A3,B3,C3,D3] = mergeParlNetw(A1,B1,C1,D1,A2,B2,C2,D2); % merge models
[A,B,C,D] = elimBus(A3,B3,C3,D3,subset); % keep only inputs/outputs connected to buses 1,2,3

% Prepare infinite-bus model parameters (bus 1)
v1_d = (2^0.5)*1e3*unit_params{1}.V; % [V] voltage source d-component
v1_q = 0; % [V] voltage source q-component
v1_0 = 0; % [V] voltage source 0-component

% Prepare synchronous generator model parameters (buses 2,3)
Ert = 1e3*[unit_params{2}.V, unit_params{3}.V]; % [Vrms] no-load output voltage at frequency w=ws
Pg = 1e6*[unit_params{2}.P, unit_params{3}.P]; % [W] nominal output power (single phase)
poles = 2; % number of machine poles (must be even)
J = 0.0125*Pg/ws; % [kg*m^2] rotor moment of inertia
Kd = 20*J*ws; % [W*s] damping factor

% Build routing vectors for connecting the unit models to the network
Nu = length(subset); % number of buses connected to units/generators/loads
units_to_net = [1:3:3*Nu , 2:3:3*Nu , 3:3:3*Nu];
[~,net_to_units] = sort(units_to_net);

%% Transient Simulation
% a) Run the Simulink model, and observe dq0 and abc signals
% b) Toggle the manual switch to create a step in the
% mechanical power on bus 2. Observe the change in 
% operating point, and the resulting transients in
% frequencies and powers.

%% Linearization & analysis
Tutorial_1_sim; % open Simulink model
io = getlinio(bdroot); % get linear model inputs/outputs
linsys = linearize(bdroot,io,3); % linearize the system at the operating point
[Ap,Bp,Cp,Dp]=ssdata(linsys); % the resulting linear system.
figure;  bode(linsys);  % Bode plot
figure;  pzmap(linsys); % pole map
