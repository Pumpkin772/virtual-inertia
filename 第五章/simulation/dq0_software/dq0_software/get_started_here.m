% This script provides a short tutorial
% demonstrating several basic functions of
% this software package.
%
% The dynamic model of a small test-case system
% is constructed, based on dq0 signals.
% The resulting model is provided as a Matlab
% state-space object.
%
% For more advanced examples, go to folder examples and 
% try 'ex_get_started_here.m' and other examples
%
% THE EXAMPLE IN THIS FILE
% The network includes three buses:
% Bus 1 is connected to bus 2 through inductance L12 and resistor R12
% Bus 3 is connected to bus 2 through inductance L23 and resistor R23
% Bus 2 has a resistive shunt element to ground, with resistance R2
% The dynamic model is first constructed with symbolic variables.
% This is done to demonstrate the structure of a typical state-space model.
% The dynamic model is then constructed again using numeric values.
%
% The dynamic model of the network is based
% on dq0 signals and is given by:
% d/dt x = Ax+Bu
% y = Cx+Du
% where:
% N is the number of buses
% u = [Vd(t); Vq(t); V0(t)] (vector 3Nx1).
% The vectors Vd, Vq, V0 are the dq0 transformation of the bus voltages.
% y = [Id(t); Iq(t); I0(t)] (vector 3Nx1).
% The vectors Id, Iq, I0 are the dq0 transformation of the bus currents.
%
% For more details on modeling using dq0 signals
% see function 'ssNetw'

clc;

% System frequency:
fs = 50;  % system frequency [Hz]
ws = 2*pi*fs; % [rad/s]

% % Define the network components and topology:
% % (for more details refer to function 'ssNetw')
% syms R12 R23 L12 L23 R2;  % define network elements as symbolic variables
% R_bus = [inf R2 inf];    % shunt resistances
% L_bus = [inf inf inf]; % shunt inductances
% C_bus = [inf inf inf];  % shunt capacitances
% Rtil_bus = [inf inf inf]; % resistors in series to shunt capacitors
% Tau = sym(zeros(3)); % transformer ratios
% Rb = [0 R12 0; 0 0 R23; 0 0 0]; % branch resistances
% Lb = [0 L12 0; 0 0 L23; 0 0 0]; % branch inductances
% 
% % Construct a symbolic dynamic model:
% [Asym, Bsym, Csym, Dsym] = ssNetwSym(R_bus, L_bus, C_bus, Rtil_bus, Rb, Lb, Tau);
% 
% % Display the resulting state-space matrices
% disp('State-space model of the network:')
% Asym
% Bsym
% Csym
% Dsym

% Define the network components and topology:
% (for more details refer to function 'ssNetw')
syms R12 L12 C1 C2;  % define network elements as symbolic variables
R_bus = [inf inf ];    % shunt resistances
L_bus = [inf inf ]; % shunt inductances
C_bus = [C1 C2 ];  % shunt capacitances
Rtil_bus = [inf inf ]; % resistors in series to shunt capacitors
Tau = sym(zeros(2)); % transformer ratios
Rb = [0 R12; 0 0]; % branch resistances
Lb = [0 L12; 0 0]; % branch inductances

% Construct a symbolic dynamic model:
[Asym, Bsym, Csym, Dsym] = ssNetwSym(R_bus, L_bus, C_bus, Rtil_bus, Rb, Lb, Tau);

% Display the resulting state-space matrices
disp('State-space model of the network:')
Asym
Bsym
Csym
Dsym



% Note on different options to construct dynamic models:
% - the function 'ssNetwSym' (used above)
%       constructs the dynamic model using symbolic variables.
%       This function may be used to generate small examples,
%       and enables theoretical testing of small test-case systems.
% - the function 'ssNetw' constructs the dynamic model
%       using numeric values. This model is more suitable for large power
%       networks, and may handle networks with thousands of buses.
% - the function 'ssNetwMatPower' constructs the dynamic
%       model of a network represented by a MatPower database.

% Construct the dynamic model again using numeric values.
% Define the network components and topology:
R12 = 0.1;
R23 = 0.1;
L12 = 0.1;
L23 = 0.2;
R2 = 10;
R_bus = [inf R2 inf];    % shunt resistances
L_bus = [inf inf inf]; % shunt inductances
C_bus = [inf inf inf];  % shunt capacitances
Rtil_bus = [inf inf inf]; % resistors in series to shunt capacitors
Tau = zeros(3); % transformer ratios
Rb = [0 R12 0; 0 0 R23; 0 0 0]; % branch resistances
Lb = [0 L12 0; 0 0 L23; 0 0 0]; % branch inductances

% Construct a (numeric) dynamic model:
[A, B, C, D] = ssNetw(R_bus, L_bus, C_bus, Rtil_bus, Rb, Lb, Tau, ws);

% The symbolic and numeric models are equivalent.
% This can be tested by the following code:
% A == double(subs(Asym))
% B == double(subs(Bsym))
% C == double(subs(Csym))
% D == double(subs(Dsym))

% Create a quasi-static model.
% (For more details, refer to function 'createQS')
% The quasi-static model is a pure gain model,
% obtained by the approximation s-->0 (low frequency approximation).
% Under this approximation the system can be
% modeled by means of time-varying phasors.
[Aqs, Bqs, Cqs, Dqs] = createQS(Asym, Bsym, Csym, Dsym);

% Create the admittance matrix (Ybus):
% This function creates the admittance matrix
% using the full dq0 state-space model as input.
Ybus = createYbus(Asym, Bsym, Csym, Dsym);
disp('Admittance matrix:')
Ybus

% To generate the numeric admittance matrix use
% Ybus = createYbus(A, B, C, D);

% Construct Matlab state-space models.
% The models represent the the network dynamics
% in the dq0 reference frame.
% Convert matrices to full (instead of sparse).
Ann = full(A);
Bnn = full(B);
Cnn = full(C);
Dnn = full(D);
sysDQ0 = ss(Ann, Bnn, Cnn, Dnn); % full dq0 state-space model
% Quasi-static model:
% 'double(subs(...)' converts the symbolic expressions to numeric ones.
Aqsnn = double(subs(Aqs)); 
Bqsnn = double(subs(Bqs));
Cqsnn = double(subs(Cqs));
Dqsnn = double(subs(Dqs));
sysQS = ss(Aqsnn, Bqsnn, Cqsnn, Dqsnn); % quasi-static state-space model

% View Bode plots comparing the dq0 and quasi-static models
% in the frequency domain.
% Note: the quasi-static model is characterized
% by a constant gain.
input_set = [1 4]; % These are the d,q components of the voltage on bus 1.
output_set = [2 4]; % These are the d,q components of the current on bus 2.
figure(1);
bode(sysDQ0(output_set,input_set),...
    sysQS(output_set,input_set));
title('Blue = DQ0 model,  Red = quasi-static model');
