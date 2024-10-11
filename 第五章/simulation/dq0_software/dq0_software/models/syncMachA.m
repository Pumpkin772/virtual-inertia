function [Aphi, Bphi, Cphi, Dphi, Gphi] = syncMachA(params, op, ws)
% SYNCMACHA Constructs a linear state-space model of a synchronous machine.
% 
% Usage:
%       [Aphi, Bphi, Cphi, Dphi, Gphi] = SYNCMACHA(params, op, ws)
% 
% where
%       params - structure containing the generator data
%       params.poles - number of poles (must be even)
%       params.J - [W*s^3] rotor moment of inertia
%       params.Kd - [W*s] damping factor
%       op - structure containing the machine's operating point
%       op.Vmag - [Vrms] voltage magnitude (RMS)
%       op.ph - [rad/s] voltage phase (typically in respect to bus 1)
%       op.P - [W] active power output (of a single phase)
%       op.Q - [VAr] reactive power output (of a single phase)
%       ws - [rad/s] nominal system frequency
% 
% Outputs:
%       Aphi, Bphi, Cphi, Dphi, Gphi - system matrices


if nargin < 3
    error('syncMachA:NotEnoughInputArguments', ...
        'Not enough input arguments.');
end

% Manipulation of constants:
% Operating point:
Id = (2^0.5)*real((op.P + 1j*op.Q)/(op.Vmag*exp(1j*op.ph))); % Id at operating point
Iq = -(2^0.5)*imag((op.P + 1j*op.Q)/(op.Vmag*exp(1j*op.ph))); % Iq at operating point
% Constants:
u1bar = Id; 
u2bar = Iq;
K1 = params.poles/(2*params.J*ws);
Ve = op.Vmag*2^0.5;
x1bar = op.ph;

% Linear state-space model:
Aphi = [0 1;
    (3/2)*K1*Ve*(sin(x1bar)*u1bar - cos(x1bar)*u2bar) -K1*params.Kd];
Bphi = [0 0 0;
    -(3/2)*K1*Ve*cos(x1bar) -(3/2)*K1*Ve*sin(x1bar) 0];
Cphi = [-Ve*sin(x1bar) 0;
    Ve*cos(x1bar) 0;
    0 0];
Dphi = zeros(3,3);
Gphi = [0; K1];

end
