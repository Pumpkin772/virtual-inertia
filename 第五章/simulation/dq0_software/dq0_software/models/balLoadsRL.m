function [A, B, C, D] = balLoadsRL(R_bus, L_bus, ws)
% BALLOADSRL Generates a dynamic model of balanced three-phase series RL loads.
%
% Usage:
%       [A, B, C, D] = BALLOADSRL(R_bus, L_bus, ws)
% 
% where
%       R_bus - vector Nx1, resistance on each bus in Ohm
%       L_bus - vector Nx1, inductance on each bus in H
%       ws - the nominal system frequency [rad/s]
%
% Outputs:
%       A, B, C, D - sparse system matrices
% 
% 
% **** Comments: ****
% The dynamic model is represented in state-space form as
% d/dt x = A*x + B*u
% y = C*x + D*u
% where:
% N is the number of buses
% u = [Vd(t); Vq(t); V0(t)] (vector 3Nx1).
% The vectors Vd, Vq, V0 are the dq0 transformation of the bus voltages.
% y = [Id(t); Iq(t); I0(t)] (vector 3Nx1).
% The vectors Id, Iq, I0 are the dq0 transformation of the injected bus currents.
%
% Input Limitations:
%       All input values may not be negative.
%       R_bus[i] cannot be infinite
%       R_bus[i] and L_bus[i] cannot be both zeros
%
% Representation of elements which do not exist:
%       If L_bus[i]=0 then the inductor does not exists (resistive load)
%       If R_bus[i] = 0 then the resistor does not exists (inductive load)
%       If a load on bus i does not exist define L_bus[i] = inf, R_bus[i] = 0
% 
% See also: ssNetw


narginchk(3,3)

% Convert shunt element inputs to column vectors
R_bus = R_bus(:);
L_bus = L_bus(:);

N = length(R_bus);

% Check input dimensions
err = 0;
if (~isequal(size(R_bus),[N,1]))
    err = 1;
end
if (~isequal(size(L_bus),[N,1]))
    err = 1;
end
if (err)
    error('balLoadsRL:IncorrectDimensions', ...
        'Input dimensions must agree.');
end

% Check positive inputs
if any(R_bus < 0)
    error('balLoadsRL:IncorrectInput', ...
        'Resistance must be nonnegative.');
end
if any(L_bus < 0)
    error('balLoadsRL:IncorrectInput', ...
        'Inductance must be nonnegative.');
end
% Check that R_bus[i] and L_bus[i] are not both zeros
if any((L_bus==0) & (R_bus==0))
    error('balLoadsRL:IncorrectInput', ...
        'Resistance and inductance cannot be both zero.');
end
% Check that R_bus[i] is finite
if any(isinf(R_bus))
    error('balLoadsRL:IncorrectInput', ...
        'Resistance cannot be infinite.');
end
    
% Convert shunt element matrices to full matrices
R_bus = full(R_bus);
L_bus = full(L_bus);
Adq = sparse(2*N,2*N);
Bdq = sparse(2*N,2*N);
Cdq = sparse(2*N,2*N);
Ddq = sparse(2*N,2*N);
for ii=1:N
    L=L_bus(ii); R=R_bus(ii);
    if (L>0)
        pos = 2*ii - 1;
        Adq(pos,pos) = -R/L;
        Adq(pos,pos+1) = ws;
        Adq(pos+1,pos) = -ws;
        Adq(pos+1,pos+1) = -R/L;
        
        Bdq(pos , ii) = 1/(2*L);
        Bdq(pos+1 , ii+N) = 1/(2*L);
        
        Cdq(ii , pos) = 2;
        Cdq(ii+N , pos+1) = 2;
    else
        % L=0
        Ddq(ii,ii) = 1/R;
        Ddq(ii+N,ii+N) = 1/R;
    end
end

A0 = sparse(N,N);
B0 = sparse(N,N);
C0 = sparse(N,N);
D0 = sparse(N,N);
for ii=1:N
    L=L_bus(ii); R=R_bus(ii);
    if (L>0)
        A0(ii,ii) = -R/L;      
        B0(ii,ii) = 1/L;
        C0(ii,ii) = 1;
    else
        % L=0
        D0(ii,ii) = 1/R;
    end
end

% Merge dq and 0 components:
A = [Adq sparse(2*N,N); ...
     sparse(N,2*N) A0];
B = [Bdq sparse(2*N,N); ...
     sparse(N,2*N) B0];
C = [Cdq sparse(2*N,N); ...
     sparse(N,2*N) C0];
D = [Ddq sparse(2*N,N); ...
     sparse(N,2*N) D0];   

% Elimination of non-active states
% for which the corresponding rows in B are zero
A(:,~any(B,2)) = [];  % remove columns
A(~any(B,2),:) = [];  % remove rows
C(:,~any(B,2)) = [];  % remove columns
B(~any(B,2),:) = [];  % remove rows

end
