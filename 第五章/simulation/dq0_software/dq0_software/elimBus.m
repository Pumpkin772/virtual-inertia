function [Ar, Br, Cr, Dr] = elimBus(A, B, C, D, subset)
% ELIMBUS Eliminates disconnected buses.
%
% Usage:
%       [Ar, Br, Cr, Dr] = elimBus(A, B, C, D, subset)
% 
% where
%       A, B, C, D - matrices of the original model
%       subset - a vector of bus indices specifying the buses
%           to be included in the new model. Voltages and currents
%           of buses in this list are NOT eliminated, and will be
%           the inputs and outputs of the new model.
%
% Outputs:
%       Ar, Br, Cr, Dr - system matrices of the new dynamic model
% 
% 
% **** Comments: ****
% The need to eliminate disconnected buses arises in several occasions:
% a) Sometimes certain buses are not connected to either a generator
% or a load. In this case the current injected into the bus is zero. 
% b) Frequently loads are modeled as shunt elements and are integrated
% into the network model. In this case load buses appear as
% disconnected buses with zero current.
% c) In many scenarios there is a need to analyze the dynamics or
% stability of a certain subset of units in the network, typically
% only the generators. In such cases elimination of the disconnected
% buses results in a simpler dynamic model in which the inputs and outputs
% relate only to the required subset of buses. This is usually done
% after integrating the loads in the network model as shunt elements.
%
% The dynamic model of the original system is
% d/dt x = A*x + B*u
% y = C*x + D*u
% where:
% N is the number of buses
% u = [Vd(t); Vq(t); V0(t)] (vector 3Nx1).
% The vectors Vd, Vq, V0 are the dq0 transformation of the bus voltages.
% y = [Id(t); Iq(t); I0(t)] (vector 3Nx1).
% The vectors Id, Iq, I0 are the dq0 transformation of the bus currents.
%
% The dynamic model after bus elimination is
% d/dt x = Ar*x + Br*u
% y = Cr*x + Dr*u
% where:
% M is the number of remaining buses (equal to the length of 'subset')
% u = [Vd(t); Vq(t); V0(t)] (vector 3Mx1).
% Vd, Vq, V0 relate to voltages of the remaining buses.
% y = [Id(t); Iq(t); I0(t)] (vector 3Mx1).
% Id, Iq, I0 relate to currents of the remaining buses.

narginchk(5,5)

if (isempty(subset))
    Ar = [];
    Br = [];
    Cr = [];
    Dr = [];
    return;
end 

N = size(B,2)/3; % number of buses in the original model
M = length(subset); % number of buses in the new model
subset = subset(:); % convert to a column vector
subset = sort(subset); % sort from low to high

input_symbolic = isa(A,'sym') && ...
    isa(B,'sym') && ...
    isa(C,'sym') && ...
    isa(D,'sym');

diagD = diag(D);

% Check that D is diagonal
if (~isempty(find(D - diag(diagD), 1)))   
    error('elimBus:InvalidStructure', ...
        'Matrix D must be diagonal.');
end

% Check that subset contains unique values
if (length(unique(subset)) ~= length(subset))
    error('elimBus:InvalidStructure', ...
        'Data in input subset is not unique. Every bus should appear only once.');
end

% Indices of entries to keep in input/output vectors:
kp_ind = [subset ; subset+N ; subset+2*N]; 
tmp = ones(3*N,1);
tmp(kp_ind) = 0;
%  Indices of entries to eliminate:
elt_ind = find(tmp);

% Zeros on the diagonal of D
ind_d0 = find(diagD==0);

% Indices of inputs/outputs to eliminate with D(i,i)=0
% and non-zero rows in matrix C (Ci=0)
tmp = intersect(elt_ind,ind_d0);
inter_d0 = intersect(find(~all(C==0,2)),tmp);

% Indices of inputs/outputs to eliminate with D(i,i)=0
% and zero rows in matrix C (Ci=0)
zzrow = intersect(find(all(C==0,2)),tmp);

% Handling zero values on the diagonal of D:
if (~isempty(inter_d0))
    if (input_symbolic)
        disp('Error in elimBus: eliminated bus has zero shunt resistance.');
        disp('This is not allowed with symbolic inputs.');
        disp('To fix this error add a shunt resistance on the bus being eliminated.');
        return;
    else             
% The function transforms the state variables and compute a new equivalent state-space model
% (with fewer states) in which the rows Ci are zeroed.
        Ctag = C(inter_d0,:);
        if (~isempty(Ctag))
% The matrix Ctag defines constraints on the state vector: Ctag*x = 0
            Ctag = sparse(Ctag); % convert to a sparse matrix (just in case...)
                         
% The following code transforms the state vector x to a new
% (reduced order) state vector named psi2, for which Ctag*x=0
            
% We first compute a LU decomposition of Ctag:
% PP*Ctag*QQ = LL*UU. A few algebric manipulations:
% Ctag*x = 0 --> PP*Ctag*x = 0
% now define a new state vector psi such that x=QQ*psi, then
% PP*Ctag*QQ*psi = 0
% using the LU decomposition above:
% LL*UU*psi = 0, and since LL is full rank, UU*psi = 0
% UU is now divided to two matrices UU1, UU2 such that UU1 is rectangular
 % [UU1 UU2]*[psi1 ; psi2] = 0
            [LL,UU,~,QQ] = lu(Ctag); % LU decomposition
            UU1 = UU(:,1:size(UU,1));
            UU2 = UU(:,(size(UU,1)+1):end);
% Test resulting matrices:
            theps = 1e6;
            if ((condest(UU1)>theps) || (condest(QQ)>theps) ...
                    || (condest(LL)>theps))
                disp('Error in elimBus: LU decomposition failed for matrix Ctag');
                disp('Resulting matrices are badly scaled or cannot be inverted');
                disp('To fix this error add a shunt resistance on the bus being eliminated.');
                return;
            end
            
% Define HH = -inv(UU1)*UU2 to obtain
% psi1 = HH*psi2
% also define matrix SS=[HH; eye] such that
% psi = [psi1 ; psi2] = SS*psi2.
            invQQ = inv(QQ);  % x = Q*psi ,  psi=inv(Q)*x
            HH = -inv(UU1)*UU2;  % psi1 = HH*psi2
            SS = [HH ; speye(size(UU2,2))];  % psi = SS*psi2
            
% Define Wa, Wb such that psi1=Wa*psi, psi2=Wb*psi.
            Wa = [speye(size(HH,1)) , sparse(size(HH,1),size(QQ,1)-size(HH,1)) ];
            Wb = speye(size(QQ,1));
            Wb(1:size(UU1,1),:)=[];   % psi2=Wb*psi
            
% These results leads to a new state space
% (d/dt)(psi2) = Anew*psi2 + Bnew*u
% y = Cnew*psi2 + D
% (u is the original input, and y is the original output)
% in which:
            Anew = Wb*invQQ*A*QQ*SS;
            Bnew = Wb*invQQ*B;
            Cnew = C*QQ*SS;
            
% At this point, if D(i,i)=0
% than the i'th row in the matrix Cnew should be zero as well.
% test that this is indeed true:
            th = 1e-10;
            if (max(max(abs(  Cnew(inter_d0,:)  ))) > th)
                disp('Error in elimBus: failed to zero rows in matrix C');
                disp('To fix this error add a shunt resistance on the bus being eliminated.');
                return;
            end
            Cnew(inter_d0,:) = 0; % remove small numeric values (this makes sure that zeros are accurate)
            
% Derive constraints on input vector.
% These constraint arise from the dynamics of psi1
% which must fulfill the condition
% (d/dt)psi1 = HH*(d/dt)psi2
% the resulting constraints is defined as T1a*u = T1b*psi2.
% These matrices can be computed from the expressions above,
% and result in
            T1a = (Wa-HH*Wb)*invQQ*B;
            T1b = -(Wa-HH*Wb)*invQQ*A*QQ*SS;
            
% Sort state variables to appear in the same order as
% in the original state vector.
% (This part has no effect on the model dynamic behaviour,
% and is done to enable straight-forward comparison with other models).
% To this end the state variables are
% re-ordered by the transformation x=MT*psi2
% (x is the new state vector)
            zop = Wb*invQQ;
            [ii,~] = find(zop);
            sizw = size(Wb,1);
            if (sizw ~= length(ii))
                disp('Error in elimBus: failed to transform variables');
                return;
            end
            [~,iis] = sort(ii);
            MT = sparse(sizw,sizw);
            MT(sub2ind([sizw,sizw],1:sizw,iis')) = 1;          
            invMT = inv(MT);
            
% The new state-space model with state vector x
            A = invMT*Anew*MT;
            B = invMT*Bnew;
            C = Cnew*MT;
% Also sort state-variables of the constraint according to the order in x
            T1b = T1b*invMT;
        end
    end
end
% Create empty input constraints if they do not exist:
if ~exist('T1a') 
    T1a = sparse(0,size(B,2));
    T1b = sparse(0,size(B,1));
end

% The objective of the following code is to use the constraints
% on the original input vector u to express it as
% original_input = u = F1*x + F2*(desired subset of inputs)
% (x is the state-vector)

% The first step toward that end is to gather all the constraints
% on the input vector in matrix form:
% K*u = [T1a ; T2a ; Gtag ; T4a]*u
% = [T1b*x, T2b*x , (desired input subset) , 0]
% T1a, T1b are constraints that originates from the transformation of
% variables above.
% T2a, T2b are constraints of the form 0=Ci*x+Di*ui (for Di~=0)
% Gtag relate original inputs to desired set of inputs, such that
% Gtag*u = (desired input subset)
% T4a are constraints on inputs for which Di=0 and the corresponding
% row in the original matrix C (before transformation of variables)
% is also zero. In this case u(i)=0.

% iIndices of inputs/outputs to eliminate with D(i,i)~=0:
inter_d0_nz = intersect(elt_ind,find(diagD~=0)); 

T2a = sparse(length(inter_d0_nz),3*N);
T2a(sub2ind(size(T2a),1:length(inter_d0_nz),inter_d0_nz')) = 1;
if (input_symbolic || isempty(inter_d0_nz))  
    T2b = -diag(1./diagD(inter_d0_nz))*C(inter_d0_nz,:);  
else
    T2b = -diag(sparse(1./diagD(inter_d0_nz)))*C(inter_d0_nz,:);
end
Gtag = speye(3*N);
Gtag(elt_ind,:) = [];
T4a = sparse(length(zzrow),3*N);
T4a(sub2ind(size(T4a),1:length(zzrow),zzrow')) = 1;
K = [T1a; T2a; Gtag; T4a];
if (condest(K) > 1e9)
    disp('Error in elimBus: matrix K is badly scaled');
    disp('To fix this error add a shunt resistance on the bus being eliminated.');
    return;
end

% Using the equations above, the original input vector u is now expressed as
% u = F1*x + F2*(desired subset of inputs)
invK = inv(K);
TT3 = [T1b; T2b];
F1 = invK(:,1:size(TT3,1))*TT3;
F2 = invK(:,(size(TT3,1)+1):(size(TT3,1)+3*M));

% The input vector is substituted in the
% state equations to derive a new dynamic model
% with eliminated buses. The new model is
% (d/dt)x = Ar*x + Br*(desired subset of inputs)
% (desired subset of outputs) = Cr*x + Dr*(desired subset of inputs)
% Substitution of the equations above, and using the relation
% (desired subset of outputs)=Gtag*y leads to
Ar = A + B*F1;
Br = B*F2;
Cr = Gtag*(C + D*F1);
Dr = Gtag*D*F2;

end
