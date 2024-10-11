function [Aksi, Bksi, Cksi, Dksi, YbusPU] = ssNetwMatPower(mpc, ws)
% SSNETWMATPOWER Generates a dynamic model of a three-phase balanced or
% symmetrically configured transmission network represented by a
% MatPower database (mpc).
% 
% Usage:
%       [Aksi, Bksi, Cksi, Dksi] = SSNETWMATPOWER(mpc, ws)
%
% where
%       mpc - Matpower database definig network topology
%       ws - nominal system frequency [rad/s]
%
% Outputs:
%       Aksi, Bksi, Cksi, Dksi - sparse system matrices
%       YbusPU - admittance matrix in per-unit, computed using
%           the dq0 model (See function 'createYbus').
%           This matrix is useful for validating the
%           resulting model, since Ybus_pu should be identical
%           to the Ybus matrix generated directly by Matpower
%           using the 'makeYbus' function. If the model is correct
%           both matrices should be nearly equal.
% 
% See also: ssNetw, ssNetwSym, createYbus
% 
% Written by Dr. Yoash Levron, Dr. Juri Belikov, June 2016


narginchk(2,2)

% Define standard column numbers used in the MatPower database
BUS_I = 1; BUS_TYPE = 2; PD = 3; QD = 4; GS = 5; BS = 6; VM = 8;
VA = 9; BASE_KV = 10; GEN_BUS = 1; PG = 2; QG = 3; VG = 6; MBASE = 7;
GEN_STATUS = 8; F_BUS = 1; T_BUS = 2; BR_R = 3; BR_X = 4; BR_B = 5;
TAP = 9; SHIFT = 10; BR_STATUS = 11; PF = 14; QF = 15; PT = 16; QT = 17;

% Per-unit base:
baseP = 1e6*mpc.baseMVA; % [W] per-unit base power.
baseV = 1e3*mpc.bus(:,BASE_KV); % [V] per-unit base voltage
baseV(baseV==0) = 1; % sometimes the base voltage is not defined.
baseP(baseP==1e9) = 1;
baseI = baseP./baseV; % [A]
baseZ = (baseV.^2)./baseP; % [ohm]

gs = mpc.bus(:,GS); % shunt conductance (MW demanded at V = 1.0 p.u.)
bs = mpc.bus(:,BS); % shunt susceptance (MVAr injected at V = 1.0 p.u.)
from = mpc.branch(:,F_BUS); % "from" bus number
to = mpc.branch(:,T_BUS); % "to" bus number
r = mpc.branch(:,BR_R); % resistance (p.u.)
x = mpc.branch(:,BR_X); % reactance (p.u.)
b = mpc.branch(:,BR_B); % total line charging susceptance (p.u.)
tap = mpc.branch(:,TAP); % transformer tap ratio
shift = mpc.branch(:,SHIFT); % transformer phase shift angle (degrees)
stat = mpc.branch(:,BR_STATUS); % 1 = in-service, 0 = out-of-service

% Convert zeros to ones in transformer taps:
tap(tap==0) = 1;

% Test input:
if (any(gs<0))
    disp('Error in ssNetwMatPower: mpc.bus(i,GS)<0 represents.');
    disp('Active power injection (a generator) and cannot be modeled by a shunt element.');
    disp('Define mpc.bus(i,GS) to be zero or positive.');
    return;
end
% Search for duplicate branches
mnft = min([from,to],[],2);
mxft = max([from,to],[],2);
indft = [mnft , mxft];
[~,~,ic] = unique(indft,'rows');
[val,ind2] = sort(ic);
ind33 = find(diff(val) == 0);
if (~isempty(ind33))
    b_indices = ind2([ind33(1), ind33(1)+1]);
    disp(sprintf('Error in statespace_dq0_matpower: \nthe network contain a duplicate branch\nin indices %d, %d.\nPlease remove the duplicate\nby typing: mpc.branch(%d,:)=[];',b_indices(1),b_indices(2),b_indices(2)));
    disp(' ');
    return;
end
% Search for phase-shifting transformers:
if (any(shift))
    disp('Warning in statespace_dq0_matpower: transformer phase-shifts are ignored in dynamic modeling.');
    disp('Define mpc.branch(:,SHIFT)=0.');
end

N = length(gs); % number of buses
brncN = length(from); % number of branches

% Initialize network data:
R_bus = sparse(N,1);
L_bus = sparse(N,1);
C_bus = sparse(N,1);
Rtil_bus = sparse(N,1);
Rb = sparse(N,N);
Lb = sparse(N,N);
Tau = sparse(N,N);

% Scan branches and construct network data:
for nv=1:brncN
    if (stat(nv)==1) % if line is in service
        mn = min(from(nv),to(nv));
        mx = max(from(nv),to(nv));
% Physical transformers appear in SI units:
        ratio = tap(nv)*baseV(from(nv))/baseV(to(nv)); 
        Lb(mn,mx) = x(nv) * baseZ(to(nv))/ws;
        Rb(mn,mx) = r(nv) * baseZ(to(nv));
        Tau(from(nv),to(nv)) = ratio;
        
        cfrom = b(nv)*(1/(2*ws))/baseZ(to(nv));
        cfrom = cfrom/(ratio^2);
        cto = b(nv)*(1/(2*ws))/baseZ(to(nv));
        C_bus(from(nv)) = C_bus(from(nv)) + cfrom;
        C_bus(to(nv)) = C_bus(to(nv)) + cto;
    end
end
% Scan bus shunt elements
for nv=1:N    
    if ((gs(nv)>=0) && (bs(nv)>=0))
% Capacitive load
        res = baseV(nv)^2/(1e6*gs(nv));
        if ( ~isinf(res) )
            R_bus(nv) = res;
        end             
        C_bus(nv) = C_bus(nv) + (1/ws)*(1e6*bs(nv))/baseV(nv)^2;
    elseif ( (gs(nv)>=0) && (bs(nv)<0) )
% Inductive load
        res = baseV(nv)^2/(1e6*gs(nv));
        if ( ~isinf(res) )
            R_bus(nv) = res;
        end
        L_bus(nv) = (1/ws)*baseV(nv)^2/(-1e6*bs(nv));
    end
end
% Define Rtil_bus:
factor = 1e5; % RC zero frequency in comparison to Ws 1e4
ii = find(C_bus);
Rtil_bus(ii) = 1./(factor*ws*C_bus(ii));

% Build a minimal state-space model
[Aksi, Bksi, Cksi, Dksi] = ssNetw(R_bus, L_bus, C_bus, Rtil_bus, Rb, Lb, Tau, ws);

% Extract Ybus matrix from state-space model:
Ybus = createYbus(Aksi,Bksi,Cksi,Dksi);
if (isnan(Ybus))
    error('ssNetwMatPower:Fail', ...
        'Ybus cannot be generated.');
end
% Convert Ybus to per-unit:
indvec = find(Ybus);
[iivec,jjvec] = ind2sub(size(Ybus),indvec);
YbusPU = sparse(size(Ybus,1),size(Ybus,2));
for pp = 1:length(indvec)
    ind = indvec(pp);
    ii = iivec(pp);
    jj = jjvec(pp);
    YbusPU(ii,jj) = Ybus(ii,jj)*(baseV(ii)*baseV(jj)/baseP);
end

end
