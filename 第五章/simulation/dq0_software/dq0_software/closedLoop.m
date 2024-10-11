function [A, B, C, D] = closedLoop(Aksi, Bksi, Cksi, Dksi, unitData)
% CLOSEDLOOP Constructs a linear state-space model of a complete power system:
% transmission network, generators, and loads.

% Usage:
%       [A, B, C, D] = closedLoop(Aksi, Bksi, Cksi, Dksi, unitData)
% 
% where
%       Aksi, Bksi, Cksi, Dksi - system matrices of the network model
%       unitData - a cell array (Nx1) that stores the system matrices of the
%           units. Each cell unitData{n} is a structure with the following fields:
%       unitData{n}.A = matrix An
%       unitData{n}.B = matrix Bn
%       unitData{n}.C = matrix Cn
%       unitData{n}.D = matrix Dn
%       unitData{n}.G = matrix Gn
%
% Outputs:
%       A, B, C, D - system matrices of the complete power system
% 
% 
% **** Comments: ****
% All dynamic models are given in the state-space form.
% Models of the network and various units can be created automatically
% by other functions in this software package.
% 
% The network model is:
% d/dt ksi = Anet*ksi + Bnet*u
% y = Cnet*ksi + Dnet*u
% where
% ksi is the state vector
% N is the number of buses in the network
% The model inputs are: u = [Vd(t); Vq(t); V0(t)] (vector 3Nx1)
% The model outputs are: y = [Id(t); Iq(t); I0(t)] (vector 3Nx1)
%
% Note: In all models currents are considered positive
% when flowing from a unit (generator/load) into the network.
%
% Model of an nth unit (generator/load/etc connected to the nth bus):
% d/dt psi_n = An*psi_n + Bn*[id; iq; i0] + Gn*w_n
% [vd; vq; v0] = Cn*psi_n + Dn*[id; iq; i0]
% where
% psi_n is the state vector
% The model inputs are: [id; iq; i0] (vector 3x1) - dq0 components of the unit currents
% w_n - external inputs. For instance, for a synchronous machine
% typically w_n = Pm (scalar representing mechanical power)
% If there are no external inputs then Gn=[] (an empty matrix)
% The model outputs are: [vd; vq; v0] (vector 3x1) - dq0 components of the unit voltages
%
% The complete power system (output of this function) is modeled by
% d/dt x = A*x + B*W
% [u; y] = C*x + D*W
% where
% x is the state vector
% Inputs are: W - the external inputs of all units.
% W is a concatenation of the individual inputs
% w_n defined above, such that W = [w_1; ...; w_N].
% The model output is a vector [u; y] of size 6Nx1, in which:
% u = [Vd(t); Vq(t); V0(t)] (vector 3Nx1) - network voltages, as defined above.
% y = [Id(t); Iq(t); I0(t)] (vector 3Nx1) - network currents, as defined above.
% 
% See also: ssNetw


narginchk(5,5)

N = size(Bksi,2)/3; % number of buses

% Create routing matrices R, invR:
ind = (reshape(1:3*N,N,3))';
ind = ind(:);
I1 = speye(3*N); % sparse unity matrix
R = I1(ind,:);
invR = I1(:,ind);  % note that invR == inv(R)
% Construct state-space model:
inp = cell(N,1);
for ii = 1:N
    inp{ii} = unitData{ii}.A;
end
Aphi = blkdiag(inp{:});

inp = cell(N,1);
for ii = 1:N
    inp{ii} = unitData{ii}.B;
end
Bphi = blkdiag(inp{:})*R;

inp = cell(N,1);
for ii = 1:N
    inp{ii} = unitData{ii}.C;
end
Cphi = invR*blkdiag(inp{:});

inp = cell(N,1);
for ii = 1:N
    inp{ii} = unitData{ii}.D;
end
Dphi = invR*blkdiag(inp{:})*R;

inp = cell(N,1);
for ii = 1:N
    inp{ii} = unitData{ii}.G;
end
G = blkdiag(inp{:});

% Construct system matrices of the complete system
curlyD = inv(speye(3*N) - Dksi*Dphi);
curlyE = inv(speye(3*N) - Dphi*Dksi);
A = [Aksi+Bksi*Dphi*(curlyD\Cksi) Bksi*Cphi+Bksi*Dphi*(curlyD\Dksi)*Cphi;
    Bphi*(curlyD\Cksi) Aphi+Bphi*(curlyD\Dksi)*Cphi];
B = [sparse(size(Aksi,1),size(G,2)); G];
C = [(curlyE\Dphi)*Cksi curlyE\Cphi;
     curlyD\Cksi (curlyD\Dksi)*Cphi];
D = sparse(size(C,1), size(G,2));

end
