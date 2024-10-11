function [Aksi, Bksi, Cksi, Dksi, str] = ssNetw(R_bus, L_bus, C_bus, Rtil_bus, Rb, Lb, Tau, ws)
% SSNETW Generates a dynamic model of a three-phase balanced or
% symmetrically configured transmission network.
%
% Usage:
%       [Aksi, Bksi, Cksi, Dksi] =
%                   ssNetw(R_bus, L_bus, C_bus, Rtil_bus, Rb, Lb, Tau, ws)
% 
% where
%       R_bus - vector Nx1, shunt resistance on each bus in Ohm
%       L_bus - vector Nx1, shunt inductance on each bus in H
%       C_bus - vector Nx1, shunt capacitance on each bus in F
%       Rtil_bus - vector Nx1, resistance in series with the shunt capacitance in Ohm
%       Rb - matrix NxN, the [i,k] element is the resistance in Ohm on brach connecting
%         bus i and bus k. The algorithm only uses elements in the (strictly) upper
%         triangle of this matrix.
%         Rb can be defined as a sparse matrix
%       Lb - matrix NxN, the [i,k] element is the inductance (in H) on brach connecting
%         bus i and bus k. The algorithm only uses elements in the (strictly) upper
%         triangle of this matrix.
%         Lb can be defined as a sparse matrix
%       Tau - matrix NxN. This matrix indicates the branch transformer ratios. If
%         Tau[i,k] is not 0 or 1 then the line i-->k contain a transformer with
%         ratio Tau[i,k]:1. Tau[i,k]==0 is identical to Tau[i,k]==1, and both
%         values represent a transformer with 1:1 ratio, or simply a short (no
%         transformer). A non-zero (and non-unity) entry Tau[i,k] indicates that a transformer
%         is connected on branch i--k. The transformer is always located on bus i, such that the
%         unity side of the transformer faces the line inductance and resistance. 
%         The main diagonal of Tau is ignored.
%         Notice that if Tau[i,k] is different than 0 or 1 then Tau[k,i] must be zero.
%         It is generally recommended to define tau as a sparse matrix.
%         If a branch connecting i and k does not exist then Tau[i,k] is not used
%         and can be set to zero. If the network contains no transformers then
%         Tau may be defined as a zero matrix (preferably a sparse one)
%       ws - nominal system frequency [rad/s]
% 
% Outputs:
%       Aksi, Bksi, Cksi, Dksi - sparse system matrices
% 
% 
% **** Comments: ****
% The constructed model is
% d/dt ksi = Aksi*ksi + Bksi*u
% y = Cksi*ksi + Dksi*u
% where:
% N is the number of buses
% u = [Vd(t); Vq(t); V0(t)] (vector 3Nx1).
% The vectors Vd, Vq, V0 are the dq0 transformations of the bus voltages.
% y = [Id(t); Iq(t); I0(t)] (vector 3Nx1).
% The vectors Id, Iq, I0 are the dq0 transformations of the bus currents.
%
% Input conversions:       
%       R_bus[i] = 0 is interpreted as R_bus[i] = inf (no resistor)
%       L_bus[i] = 0 is interpreted as L_bus[i] = inf (no inductor)
%       C_bus[i] = 0 is interpreted as C_bus[i] = inf(capacitor is short-circuit)
%       Rtil_bus[i] = 0 is interpreted as Rtil_bus[i] = inf (no resistor)
%       Rb[i,k] = 0 is NOT converted to infinity
%       Lb[i,k] = 0 is interpreted as Lb[i,k] = inf (no inductor, branch is disconnected)
%       Tau[i,k] = 0 is interpreted as Tau[i,k] = 1 (no transformer)
%
% Input Limitations:
%       All values except Tau[i,k] cannot be negative
%       Rb[i,k] cannot be infinite
%       If Tau[i,k] is different than 0 or 1, then Tau[k,i] must be zero
%
% Representation of elements which do not exist:
% if a branch [i,k] does not exist define
%       Rb[i,k] = 0, Lb[i,k] = 0, Tau[i,k] = 0, Tau[k,i] = 0
%       Meaning: infinite inductance, zero resistance, and transformer
%       ratio of unity, which are equivalent to an open-circuit
% if a shunt capacitor-resistor does not exist define
%       C_bus[i] = 0, Rtil_bus[i] = 0
%       Meaning: infinite capacitance and resistance, equivalent to an open-circuit
% if a shunt resistor does not exist define
%       R_bus[i] = 0
%       Meaning: infinite resistance, equivalent to an open-circuit
% if a shunt inductor does not exist define
%       L_bus[i] = 0
%       Meaning: infinite inductance, equivalent to an open-circuit
% If a transformer on branch [i,k] does not exists define
%       Tau[i,k] = 0, Tau[k,i] = 0
%       Meaning: transformer ratio of unity,equivalent to a short-circuit
% 
% See also: ssNetwSym, ssNetwMatPower
% 
% Written by Dr. Yoash Levron, Dr. Juri Belikov, June 2016


narginchk(8,8)

% Convert shunt element inputs to column vectors
R_bus = R_bus(:);
L_bus = L_bus(:);
C_bus = C_bus(:);
Rtil_bus = Rtil_bus(:);
N = length(C_bus);

% Check input dimensions
err = 0;
if (~isequal(size(R_bus),[N,1]))
    err = 1;
end
if (~isequal(size(Rtil_bus),[N,1]))
    err = 1;
end
if (~isequal(size(L_bus),[N,1]))
    err = 1;
end
if (~isequal(size(Rb),[N, N]))
    err = 1;
end
if (~isequal(size(Lb),[N, N]))
    err = 1;
end
if (~isequal(size(Tau),[N, N]))
    err = 1;
end
if (~isequal(size(ws),[1 1]))
    err = 1;
end
if (err)
    error('ssNetw:InputDimensions', ...
        'Wrong input dimensions.');
end
    
% Convert shunt element matrices to full matrices
R_bus = full(R_bus);
L_bus = full(L_bus);
C_bus = full(C_bus);
Rtil_bus = full(Rtil_bus);

% Convert zeros to infs in shunt elements:
R_bus(R_bus == 0) = inf;
L_bus(L_bus == 0) = inf;
C_bus(C_bus == 0) = inf;
Rtil_bus(Rtil_bus == 0) = inf;

% Zero the lower diagonals of Lb, Rb
Rb = triu(Rb,1);
Lb = triu(Lb,1);

% Convert infs to zeros in Lb
Lb(isinf(Lb)) = 0;

% Replace ones by zeros in matrix Tau
Tau(Tau == 1) = 0;
% Remove the main diagonal in Tau
Tau(sub2ind(size(Tau),1:N,1:N)) = 0;

% Zero transpose elements in Tau
[ii, kk] = find(triu(Tau));
Tau(sub2ind(size(Tau),kk,ii)) = 0;
[ii, kk] = find(tril(Tau));
Tau(sub2ind(size(Tau),kk,ii)) = 0;

% Create the matrix At
Ap_vec = zeros(N*(N+3),1);
ii = 1:N;
Ap_vec(2*ii-1) = -1./(C_bus.*Rtil_bus)-1j*ws;
Ap_vec(2*ii) = -1./(C_bus.*Rtil_bus)+1j*ws;
Ap_vec(2*ii-1+2*N) = -1j*ws;
Ap_vec(2*ii+2*N) = +1j*ws;
vec1 = zeros((N-1)*N/2,1);
indpp = find(Lb>0);
if (~isempty(indpp))
    [ii,kk] = ind2sub(size(Lb),indpp);
    cnt = (kk-ii)+(ii-1).*(N-ii/2);
    vec1(cnt) = -Rb(indpp)./Lb(indpp);
end
Ap_vec((4*N+1):2:end) = vec1 - 1j*ws;
Ap_vec((4*N+2):2:end) = vec1 + 1j*ws;
At = spdiags(Ap_vec,0,N*(N+3),N*(N+3));

clear vec1 Ap_vec

BL = sparse(2*N,2*N);
for ii=1:N
    ind1 = 2*ii-1;
    BL(ind1,ii) = -1/(2*C_bus(ii).*(Rtil_bus(ii).^2));
    BL(ind1+1,ii) = -1/(2*C_bus(ii).*(Rtil_bus(ii).^2));
    BL(ind1,ii+N) = -1j/(2*C_bus(ii).*(Rtil_bus(ii).^2));
    BL(ind1+1,ii+N) = +1j/(2*C_bus(ii).*(Rtil_bus(ii).^2));    
end

BR = sparse(N*(N-1),2*N);
indpp = find(Lb>0);
if (~isempty(indpp))
    [i1,i2] = ind2sub(size(Lb),indpp);
    cnt = 2*(  (i2-i1)+(i1-1).*(N-i1/2)  ) -1;
    BR(sub2ind(size(BR),cnt,i1)) = 1./(2*Lb(indpp));
    BR(sub2ind(size(BR),cnt,i2)) = -1./(2*Lb(indpp));
    BR(sub2ind(size(BR),cnt,i1+N)) = 1j./(2*Lb(indpp));
    BR(sub2ind(size(BR),cnt,i2+N)) = -1j./(2*Lb(indpp));
    BR(sub2ind(size(BR),cnt+1,i1)) = 1./(2*Lb(indpp));
    BR(sub2ind(size(BR),cnt+1,i2)) = -1./(2*Lb(indpp));
    BR(sub2ind(size(BR),cnt+1,i1+N)) = -1j./(2*Lb(indpp));
    BR(sub2ind(size(BR),cnt+1,i2+N)) = 1j./(2*Lb(indpp));
end

% The following code handles transformer ratios (Tau) in matrix BR
indpp = find(triu(Tau,1) ~= 0);    % search in the upper diagonal
if (~isempty(indpp))
    [i1,i2] = ind2sub(size(Tau),indpp);
    cnt = 2*((i2-i1)+(i1-1).*(N-i1/2))-1;
    BR(sub2ind(size(BR),cnt,i1)) = BR(sub2ind(size(BR),cnt,i1)) ./ Tau(indpp);
    BR(sub2ind(size(BR),cnt,i1+N)) = BR(sub2ind(size(BR),cnt,i1+N)) ./ Tau(indpp);
    BR(sub2ind(size(BR),cnt+1,i1)) =  BR(sub2ind(size(BR),cnt+1,i1)) ./ Tau(indpp);
    BR(sub2ind(size(BR),cnt+1,i1+N)) =  BR(sub2ind(size(BR),cnt+1,i1+N)) ./ Tau(indpp);
end
indpp = find( tril(Tau,-1) ~= 0 );    % search in the lower diagonal
if (~isempty(indpp))
    [i1,i2] = ind2sub(size(Tau),indpp);
    cnt = 2*(  (i1-i2)+(i2-1).*(N-i2/2)  ) -1;
    BR(sub2ind(size(BR),cnt,i1)) = BR(sub2ind(size(BR),cnt,i1)) ./ Tau(indpp);
    BR(sub2ind(size(BR),cnt,i1+N)) = BR(sub2ind(size(BR),cnt,i1+N)) ./ Tau(indpp);
    BR(sub2ind(size(BR),cnt+1,i1)) =  BR(sub2ind(size(BR),cnt+1,i1)) ./ Tau(indpp);
    BR(sub2ind(size(BR),cnt+1,i1+N)) =  BR(sub2ind(size(BR),cnt+1,i1+N)) ./ Tau(indpp);
end
% Done with transformers

BM = sparse(2*N,2*N);
for ii=1:N
    ind1 = 2*ii-1;
    BM(ind1,ii) = 1/(2*L_bus(ii));
    BM(ind1+1,ii) = 1/(2*L_bus(ii));
    BM(ind1,ii+N) = 1j/(2*L_bus(ii));
    BM(ind1+1,ii+N) = -1j/(2*L_bus(ii));    
end

Bt = [BL ; BM; BR];

clear BL BR BM

Ca = speye(N);
Cb = sparse(N,2*N);
Cc = sparse(N,2*N);
Cb(:,1:2:end) = Ca;
Cb(:,2:2:end) = Ca;
Cc(:,1:2:end) = (-1j)*Ca;
Cc(:,2:2:end) = 1j*Ca;
CL = [Cb ; Cc ];

clear Ca Cb Cc

CY = sparse(N,N*(N-1)/2);
ind1 = 1;  % column
ind2 = 1; % row
for jj = N-1:(-1):1
    CY(ind2,ind1:(ind1+jj-1))=1;
    CY((ind2+1):end,ind1:(ind1+jj-1)) = -speye(jj);
    ind1 = ind1+jj;
    ind2 = ind2+1;
end

CYU = sparse(N,N*(N-1));
CYD = sparse(N,N*(N-1));
if (~isempty(CY))
    CYU(:,1:2:end) = CY;
    CYU(:,2:2:end) = CY;
    CYD(:,1:2:end) = (-1j)*CY;
    CYD(:,2:2:end) = 1j*CY;
end

CR = [CYU ; CYD];
clear CYU CYD

% The following code handles transformer ratios (Tau) in matrix CR
indpp = find( triu(Tau,1) ~= 0 );    % search in the upper diagonal
if (~isempty(indpp))
    [i1,i2] = ind2sub(size(Tau),indpp);
    cnt = 2*(  (i2-i1)+(i1-1).*(N-i1/2)  ) -1;
    CR(sub2ind(size(CR),i1,cnt)) = CR(sub2ind(size(CR),i1,cnt)) ./ Tau(indpp);
    CR(sub2ind(size(CR),i1+N,cnt)) = CR(sub2ind(size(CR),i1+N,cnt)) ./ Tau(indpp);
    CR(sub2ind(size(CR),i1,cnt+1)) =  CR(sub2ind(size(CR),i1,cnt+1)) ./ Tau(indpp);
    CR(sub2ind(size(CR),i1+N,cnt+1)) =  CR(sub2ind(size(CR),i1+N,cnt+1)) ./ Tau(indpp);
end
indpp = find( tril(Tau,-1) ~= 0 );    % search in the lower diagonal
if (~isempty(indpp))
    [i1,i2] = ind2sub(size(Tau),indpp);
    cnt = 2*(  (i1-i2)+(i2-1).*(N-i2/2)  ) -1;
    CR(sub2ind(size(CR),i1,cnt)) = CR(sub2ind(size(CR),i1,cnt)) ./ Tau(indpp);
    CR(sub2ind(size(CR),i1+N,cnt)) = CR(sub2ind(size(CR),i1+N,cnt)) ./ Tau(indpp);
    CR(sub2ind(size(CR),i1,cnt+1)) =  CR(sub2ind(size(CR),i1,cnt+1)) ./ Tau(indpp);
    CR(sub2ind(size(CR),i1+N,cnt+1)) =  CR(sub2ind(size(CR),i1+N,cnt+1)) ./ Tau(indpp);
end
% Done with transformers

Ct = [CL CL CR];

clear CL CR

Dt_vec = 1./R_bus+1./Rtil_bus;
Dt_vec = [Dt_vec ; Dt_vec ; Dt_vec];
Dksi = spdiags(Dt_vec,0,3*N,3*N);

tvec1 = ones(N*(N+3),1);
tvec1(2:2:end) = -1j;
Ts = spdiags(tvec1,0,N*(N+3),N*(N+3));
tvec2 = zeros(N*(N+3)-1,1);
tvec2(1:2:end) = 1j;
Ts = Ts + [spdiags(tvec2,1,N*(N+3)-1,N*(N+3)); sparse(1,N*(N+3))];
tvec3 = zeros(N*(N+3)-1,1);
tvec3(1:2:end) = 1;
Ts = Ts + [spdiags(tvec3,-1,N*(N+3),N*(N+3)-1),  sparse(N*(N+3),1)];

tvec1 = (1/2)*ones(N*(N+3),1);
tvec1(2:2:end) = 1j/2;
invTs = spdiags(tvec1,0,N*(N+3),N*(N+3));
tvec2 = zeros(N*(N+3)-1,1);
tvec2(1:2:end) = 1/2;
invTs = invTs + [spdiags(tvec2,1,N*(N+3)-1,N*(N+3)); sparse(1,N*(N+3))];
tvec3 = zeros(N*(N+3)-1,1);
tvec3(1:2:end) = -1j/2;
invTs = invTs + [spdiags(tvec3,-1,N*(N+3),N*(N+3)-1),  sparse(N*(N+3),1)];

Adq = invTs*At*Ts;
Bdq = invTs*Bt;
Cdq = Ct*Ts;

Adq = real(Adq);
Bdq = real(Bdq);
Cdq = real(Cdq);

A0 = sparse(N*(N+3)/2,N*(N+3)/2);
for ii=1:N
    A0(ii,ii) = -1/(C_bus(ii)*Rtil_bus(ii));
end
indpp = find(Lb>0);
if (~isempty(indpp))
    [ii,kk] = ind2sub(size(Lb),indpp);
    cnt = (kk-ii)+(ii-1).*(N-ii/2);
    A0(sub2ind(size(A0),cnt+2*N,cnt+2*N)) = -Rb(indpp)./Lb(indpp);
end

B0L = sparse(N,N);
for ii=1:N
    B0L(ii,ii) = -1/(C_bus(ii)*(Rtil_bus(ii)^2));
end
B0R = sparse(N*(N-1)/2,N);
indpp = find(Lb>0);
if (~isempty(indpp))
    [i1,i2] = ind2sub(size(Lb),indpp);
    cnt = (i2-i1)+(i1-1).*(N-i1/2);
    B0R(sub2ind(size(B0R),cnt,i1)) = 1./(Lb(indpp));
    B0R(sub2ind(size(B0R),cnt,i2)) = -1./(Lb(indpp));
end

% The following code handles transformer ratios (Tau) in matrix B0R
indpp = find( triu(Tau,1) ~= 0 );    % search in the upper diagonal
if (~isempty(indpp))
    [i1,i2] = ind2sub(size(Tau),indpp);
    cnt = (i2-i1)+(i1-1).*(N-i1/2);
    B0R(sub2ind(size(B0R),cnt,i1)) = B0R(sub2ind(size(B0R),cnt,i1)) ./ Tau(indpp);
end
indpp = find( tril(Tau,-1) ~= 0 );    % search in the lower diagonal
if (~isempty(indpp))
    [i1,i2] = ind2sub(size(Tau),indpp);
    cnt = (i1-i2)+(i2-1).*(N-i2/2);
    B0R(sub2ind(size(B0R),cnt,i1)) = B0R(sub2ind(size(B0R),cnt,i1)) ./ Tau(indpp);
end
% Done with transformers

B0M = sparse(N,N);
for ii=1:N
    B0M(ii,ii) = +1/L_bus(ii);
end
B0 = [B0L ; B0M ; B0R];

CR0 = CY;

% The following code handles transformer ratios (Tau) in matrix CR0
indpp = find( triu(Tau,1) ~= 0 );    % search in the upper diagonal
if (~isempty(indpp))
    [i1,i2] = ind2sub(size(Tau),indpp);
    cnt = (i2-i1)+(i1-1).*(N-i1/2);
    CR0(sub2ind(size(CR0),i1,cnt)) = CR0(sub2ind(size(CR0),i1,cnt)) ./ Tau(indpp);
end
indpp = find( tril(Tau,-1) ~= 0 );    % search in the lower diagonal
if (~isempty(indpp))
    [i1,i2] = ind2sub(size(Tau),indpp);
    cnt = (i1-i2)+(i2-1).*(N-i2/2);
    CR0(sub2ind(size(CR0),i1,cnt)) = CR0(sub2ind(size(CR0),i1,cnt)) ./ Tau(indpp);
end
% Done with transformers

C0 = [speye(N), speye(N), CR0];

% Final Aksi, Bksi, Cksi matrices
Aksi = [Adq sparse(N*(N+3),N*(N+3)/2); ...
     sparse(N*(N+3)/2,N*(N+3)) A0];
Bksi = [Bdq sparse(N*(N+3),N); ...
     sparse(N*(N+3)/2,2*N) B0];
Cksi = [Cdq sparse(2*N,N*(N+3)/2); ...
     sparse(N,N*(N+3)), C0];

% Eliminate entries corresponding to inactive states
% zero all elements associated with redundant states.
Aksi(:, ~any(Bksi,2) ) = [];  % remove columns
Cksi(:, ~any(Bksi,2) ) = [];  % remove columns
Aksi( ~any(Bksi,2), : ) = [];  % remove rows
Bksi( ~any(Bksi,2), : ) = [];  % remove rows

end
