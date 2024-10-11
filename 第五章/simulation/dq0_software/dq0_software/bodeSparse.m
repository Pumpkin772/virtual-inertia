function [Mag, Ph] = bodeSparse(A, B, C, D, u, w)
% BODESPARSE Bode frequency response of a sparse state-space model.
%
% Usage:
%       [Mag, Ph] = bodeSparse(A, B, C, D, input, w)
%
% where
%       A, B, C, D - system matrices
%       u - input's serial number. For single input input=1.
%           For multi-input models select 1 <= input <= size(B,2)
%       w - frequency range [rad/s]
% 
% Outputs:
%       Mag - magnitude [dB]
%       Ph - phase [deg]


narginchk(6,6)
if (~isscalar(u))
    error('bodeSparse:NonScalarInput', ...
        'Input variable must be scalar.');
elseif ((u<1) || (u>size(B,2)))
    error('bodeSparse:WrongInputRange', ...
        'Input variable must be in the range 1 <= input <= size(B,2).');
end

B = B(:,u); % select input in B
D = D(:,u); % select input in D
NN = size(A,1);
speyeNN = speye(NN);
lenW = length(w);

% Compute frequency response
valvec = zeros(size(C,1),lenW);
for k=1:lenW;
    valvec(:,k) = C*((1j*w(k)*speyeNN-A)\B)+D;
end

% Compute magnitude and phase
Mag = abs(valvec);
Ph = zeros(size(C,1),lenW);
for ii=1:size(C,1)
    Ph(ii,:) = phase(valvec(ii,:))*180/pi;
end

end
