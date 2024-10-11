% Small-signal modeling and stability analysis of large power systems
% using dq0 signals.
% This file demonstrates:
% - How to model the small-signal dynamics of large power systems
%   using dq0 signals.
% - How to evaluate the system stability based on eigenvalues.
% - How to compute step responses of the complete system.
% Power systems ranging in size from 4 to 2383 buses are considered.
%
% To begin, choose an example below, and run the script.
% Note: the MatPower package should be installed (available online).
%
% This file loads test-case networks from the Matpower database
% and constructs a small-signal dynamic model of the complete
% power system (network + loads + generators) in the dq0
% reference frame. Stability at an operating point is analyzed
% by examining eigenvalues of the state matrix A,
% or by examining the step response.
%
% Model construction:
% The dynamic model is constructed in several stages:
% a) First a dynamic model of the transmission network is constructed.
% b) Loads are represented as series RL shunt elements,
%   and embedded into the network model.
% c) All buses are eliminated, expect from generator buses.
%   Following this stage a generator is connected at each bus,
%   and the number of buses is equal to the number of generators.
% d) Small-signal models of the generators are constructed.
%   The reference bus (typically bus 1) is modeled as
%   an infinite bus (a constant voltage source).
%   Other generators are modeled as synchronous machines.
%   (for details on this model please refer to 'syncMachA').
%   The models are linearized in the neighborhood of the system's
%   operating point (solution of power flow equations).
% e) The models obtained in the previous stages are combined
%   to create a unified model of the complete power system
%   (network + loads + generators).
%
% At this stage, the power system is modeled by
% d/dt x = A*x + B*W
% [u; y] = C*x + D*W
% where:
% x is the unified state vector
% W - vector of external inputs. In this case W
% is the mechanical power inputs of synchronous machines,
% and is measured in Watt. The length of W is determined by
% the number of synchronous machines in the system.
% The output vector is [u; y], where
% u = [Vd(t); Vq(t); V0(t)] vector 3Nx1 of dq0 bus voltages [V]
% y = [Id(t); Iq(t); I0(t)] vector 3Nx1 dq0 bus currents [A]
% (injected from generators to buses)
%
% Note: Following step (c) above, the number of buses N
% is equal to the number of generators.

clc;
clear all;
define_constants; % MatPower function. Defines MatPower constants.

% Choose an example from the list below:
% Stable systems examples: 1, 2, 3, 4, 6, 7, 8, 9, 10
% Unstable systems examples:  5, 11
%
% Note: example 11 demonstrates a 2383-bus network,
% and takes a few minutes to compute.

choose_example = 4;

% Load a network from the MatPower database:
switch choose_example
    case 1
        filename = 'case9'; % Matpower test-case network
        mpc = loadcase(filename); % Matpower function
    case 2
        filename = 'case_ieee30';
        mpc = loadcase(filename);
        step_Tf = 0.4; % [s] final time for step response
        step_value = 90e6; % [W] step size
    case 3
        filename = 'case57';
        mpc = loadcase(filename);
        % remove duplicate branches in the network
        mpc.branch(20,:)=[];
        mpc.branch(35,:)=[];
        step_Tf = 0.8; % [s] final time for step response
        step_value = 90e6; % [W] step size
    case 4
        % 118-bus network with stable parameters
        filename = 'case118';
        mpc = loadcase(filename);
        % remove duplicate branches
        mpc.branch(67,:)=[];
        mpc.branch(75,:)=[];
        mpc.branch(97,:)=[];
        mpc.branch(84,:)=[];
        mpc.branch(120,:)=[];
        mpc.branch(134,:)=[];
        mpc.branch(136,:)=[];
        step_Tf = 1; % [s] final time for step response
        step_value = 270e6; % [W] step size
        % synchronous machine constants:
        KD_GAIN_CNST = 200; % damping factor. Reduce this value to obtain an unstable system
        J_GAIN_CNST = 0.004; % rotor moment of inertia    
    case 5
        % 118-bus network with unstable parameters
        filename = 'case118';
        mpc = loadcase(filename);
        % remove duplicate branches
        mpc.branch(67,:)=[];
        mpc.branch(75,:)=[];
        mpc.branch(97,:)=[];
        mpc.branch(84,:)=[];
        mpc.branch(120,:)=[];
        mpc.branch(134,:)=[];
        mpc.branch(136,:)=[];
        % synchronous machine constants:
        KD_GAIN_CNST = 100; % damping factor, resulting in an unstable network
        J_GAIN_CNST = 0.004;  % rotor moment of inertia     
    case 6
        filename = 'case4gs';
        mpc = loadcase(filename);
        step_Tf = 0.8; % [s] final time for step response
        step_value = 90e6; % [W] step size.
    case 7
        filename = 'case6ww';
        mpc = loadcase(filename);
        step_Tf = 0.7; % [s] final time for step response
    case 8
        filename = 'case9Q';
        mpc = loadcase(filename);
        step_Tf = 4; % [s] final time for step response
    case 9
        filename = 'case30Q';
        mpc = loadcase(filename);
        step_Tf = 0.3; % [s] final time for step response
    case 10
        filename = 'case39';
        mpc = loadcase(filename);
        % remove loads with negative reactive power
        mpc.bus(9,QD) = 0;
        mpc.bus(24,QD) = 0;
        % determine minimal resistance of all branches (for stability)
        miniminal_branch_resistance_pu = 0.001;
        step_Tf = 0.5; % [s] final time for step response
        step_value = 300e6; % [W] step size
    case 11
        filename = 'case2383wp';  % 2383 bus network
        % Notice: this example may take several minutes to compute
        mpc = loadcase(filename); 
        % remove duplicate branches
        mpc.branch(19,:)=[];
        mpc.branch(661,:)=[];
        mpc.branch(1004,:)=[];
        mpc.branch(1675,:)=[];
        mpc.branch(2350,:)=[];
        mpc.branch(2592,:)=[];
        mpc.branch(2634,:)=[];
        mpc.branch(2789,:)=[];
        mpc.branch(2872,:)=[];
        mpc.branch(2879,:)=[];
        % remove phase shifts in transformers:
        mpc.branch(:,SHIFT)=0;
        % zero active and reactive powers of loads when negative
        ii = find(mpc.bus(:,QD) < 0);
        mpc.bus(ii,QD) = 0;
        ii = find(mpc.bus(:,PD) < 0);
        mpc.bus(ii,PD) = 0;
        miniminal_branch_resistance_pu = 5e-6;
    otherwise
        disp('Unknown Matpower test-case');
        return
end

% Process and prepare network data:
% If the per-unit base voltage is not specified, choose a default value
ind = find(mpc.bus(:,BASE_KV)==0);
mpc.bus(ind,BASE_KV)=100; % [kV]
% Sort bus data according to bus index
[~,ind1] = sort(mpc.bus(:,BUS_I));
mpc.bus = mpc.bus(ind1,:);
% Sort generator data according to generator bus index
[~,ind1] = sort(mpc.gen(:,GEN_BUS));
mpc.gen = mpc.gen(ind1,:);
% Place minimal resistances on network branches (important for stability)
if ~exist('miniminal_branch_resistance_pu') 
    miniminal_branch_resistance_pu = 0.01;
end
mpc.branch(:,BR_R) = max( mpc.branch(:,BR_R), miniminal_branch_resistance_pu ); % branch resistance p.u.

% Compute the system operating point by solving the power flow equations.
% ('runpf' is a MatPower function)
mpc = runpf(mpc);
if (mpc.success ~= 1)
    fprintf('%s\n\n',filename); % display network name
    disp('power flow failed to converge.');
    return;
end
fprintf('%s\n\n',filename); % display network name

% Define per-unit (p.u.) quantities and other constants
baseVA = mpc.baseMVA * 1e6;  % [W] 1p.u. base power
baseV = 1e3*mpc.bus(:,BASE_KV); % [V] 1p.u. base voltage
baseA = baseVA./baseV; % [A] 1p.u. base current
baseZ = baseV./baseA; % [Ohm] 1p.u. base impedance
ws = 2*pi*50; % [rad/s] system frequency
N = size(mpc.bus,1);  % number of buses in the network
fprintf('Number of buses = %d\n',N);
if (N>300)
    disp('This example may take several minutes to compute...');
    fprintf('Computing dq0 dynamic model. Please wait ...\n');
end
fprintf('\n');

% Build state-space model of the transmission network
[A1, B1, C1, D1, YbusPU] = ssNetwMatPower(mpc, ws);

% Find Bus indices of generators:
refbus = find(mpc.bus(:,BUS_TYPE) == 3); % reference bus (typically bus 1)
if (length(refbus)>1)
    disp('Error - more than one reference bus');
    return;
end
ssbus = mpc.gen(:,GEN_BUS);  % generator bus indices (including ref bus)

% Extract load parameters from the MatPower database:
% Loads are represented as series RL shunt elements
% For details refer to function 'balLoadsRL'
vv = baseV.*mpc.bus(:,VM); % [V]
pp = 1e6*mpc.bus(:,PD);  qq = 1e6*mpc.bus(:,QD); % [W]
ss2 = pp.^2 + qq.^2; % [VA]
if (any(qq<0))
    disp('Warning: loads with negative reactive power are not supported');
    return;
end
if (any(pp<0))
    disp('Warning: loads with negative active power are not supported');
    return;
end
R_bus = pp.*(vv.^2) ./ ss2;  % [Ohm] vector of load resistors
L_bus = (1/ws)*qq.*(vv.^2) ./ ss2; % [H] vector of load inductors
L_bus(ss2 == 0) = inf;   R_bus(ss2 == 0) = 0; % remove non-existing loads
% Build state-space model of loads:
[A2, B2, C2, D2] = balLoadsRL(R_bus, L_bus, ws);
% Merge state-space models of original network and load network:
[A3, B3, C3, D3] = mergeParlNetw(A1, B1, C1, D1, A2, B2, C2, D2);
% Eliminate all buses except from generator buses:
[Anet, Bnet, Cnet, Dnet] = elimBus(A3, B3, C3, D3, ssbus);

% clear A2 B2 C2 D2;
% clear A3 B3 C3 D3;

% Construct small-signal dynamic models of generators,
% and prepare the structure 'unitData' to be used
% in the function 'closedLoop'
unitData = struct([]); % create an empty structure
for nj = 1:length(ssbus)
% Construct the dynamic model of each generator.
% More details on modeling units is provided in
% functions 'closedLoop' and 'syncMachA'
    
    n = ssbus(nj); % generator bus index
% 'nj' is the generator index in mpc.gen,
% 'nj' is also the bus index in the network after bus elimination
% (the network after bus elimination includes only generator buses)
% 'n' is the bus index in the original network.
     
    if (n == refbus)
% Infinite bus (reference bus)
        [An, Bn, Cn, Dn, Gn] = infBus();
    else
% Generator is a synchronous machine
% Define synchronous machine parameters:       
        op = struct([]);
% Operating point:
        op(1).Vmag = mpc.bus(n,VM) * baseV(n); % [Vrms] voltage magnitude (RMS)
        op.ph = mpc.bus(n,VA) * pi/180; % [rad/s] voltage phase
        op.P = mpc.gen(nj,PG) * 1e6; % [W] generator active power output (of a single phase)
        op.Q = mpc.gen(nj,QG) * 1e6; % [VAr] generator reactive power output (of a single phase)
        params = struct([]);
        params(1).poles = 2; % number of machine poles
        if ~exist('KD_GAIN_CNST')
            KD_GAIN_CNST = 100;
        end
        if ~exist('J_GAIN_CNST')
            J_GAIN_CNST = 0.002;
        end      
% Rotor moment of inertia:
        Pratio = 0.2;
        if (op.P > Pratio*baseVA)
            params.J = J_GAIN_CNST*12*op.P/ws; % proportional to power
        else
            params.J = J_GAIN_CNST*12*Pratio*baseVA/ws; % not proportional to power
        end
% Machine's damping factor:
        params.Kd = KD_GAIN_CNST*params.J*ws; % damping factor
% Build synchronous machine model:
        [An, Bn, Cn, Dn, Gn] = syncMachA(params, op, ws);
    end
% add dynamic model of unit to 'unitData' structure:
    unitData{nj}.A=An; unitData{nj}.B=Bn; unitData{nj}.G=Gn;
    unitData{nj}.C=Cn; unitData{nj}.D=Dn;
end

% Construct the dynamic model of the complete power system
% (network + loads + generators) in the dq0 reference frame:
[A, B, C, D] = closedLoop(Anet, Bnet, Cnet, Dnet, unitData);
figure(3); spy(A); title('state matrix A');

% Test stability of the power system by
% examining the eigenvalues of the state matrix A
% The system is stable if and only if all eigenvalues
% have negative real parts
if (size(A,1) < 1500)
% For a relatively small matrix A, compute the eigenvalues directly
    fullA = full(A);
    eigenvalues = eig(fullA);
    clear fullA;
% Find eigenvalue with maximal real part:
    dr_max = eigenvalues( real(eigenvalues) == max(real(eigenvalues)) );
    dr_max = dr_max(1);
% Plot eigenvalues:
    figure(1);
    plot(real(eigenvalues),imag(eigenvalues),'kx','Color',[0.25 0.25 0.25],'MarkerSize',4);
    xlim([-200 max(real(dr_max)*1.1,10)]);
    title('Poles of the closed loop system');
    xlabel('Real'); ylabel('Imag');
else
% If the matrix A is very large, use 'eigs' to compute the eigenvalues.
    close all;
    number_of_tries = 20;
% Randomize several starting points and compute eigenvalues:
    dr_max=-inf;
    for ii=1:number_of_tries
        fprintf('computing eigenvalues %d/%d\n',ii,number_of_tries);
        sp = 150+100*rand(1) + 100i*randn(1);
        eg = eigs(A,2,sp);
        if (max(real(eg))>real(dr_max))
            ind = find(max(real(eg)) == real(eg));
            dr_max = eg(ind);
        end
    end
end
fprintf('Largest eigenvalue real part is %g\n',real(dr_max));
if (real(dr_max) < 0)
    fprintf('System is stable at operating point.\n\n');
else
    fprintf('System is UNSTABLE at operating point.\n\n');
    close all;
    disp('done');
    return;
end

% Plot step response of the closed loop system
input = 1; % choose system input (generator index, not including reference bus)
if (~exist('step_value'))
    step_value = 30e6; % [W] step size. (mechanical power of synchronous generator)
end
if (~exist('step_Tf'))
    step_Tf = 5; % [s] final time for step response
end
step_T = 1e-4; % [s] numeric step time for step response (sampling time)
display=0;
if (N>300)
   fprintf('Computing step response. Please wait ...\n');
   display=1;
   step_T = 1e-3;
end

% Compute and plot the step response
% The function 'stepSparse' below performs the same task as Matlab's
% function 'step', but is designed to work fast with large dynamic
% systems that are represnted by sparse system matrices. 
[y, t] = stepSparse(A, B(:,input), C, D(:,input), step_Tf, step_T, display);
y = step_value * y;
% Extract dq0 components of bus voltages and currents
% (these represent the small signal response to the step
% in the neighbourhood of the operating point)
MM = length(ssbus);
vd = y(1:MM,:);
vq = y((MM+1):2*MM,:);
v0 = y((2*MM+1):3*MM,:);
id = y((3*MM+1):4*MM,:);
iq = y((4*MM+1):5*MM,:);
i0 = y((5*MM+1):6*MM,:);
% Compute operating point in terms of dq0 components        
Vmagtt = mpc.bus(ssbus,VM).*baseV(ssbus);
phtt = mpc.bus(ssbus,VA)*pi/180;
Ptt = mpc.gen(:,PG) * 1e6; % [W] generator active power output (of a single phase)
Qtt = mpc.gen(:,QG) * 1e6; % [VAr] generator reactive power output (of a single phase)  
vd_op = (2^0.5)*real(  Vmagtt.*exp(1j*phtt)  );
vq_op = (2^0.5)*imag(  Vmagtt.*exp(1j*phtt)  );
v0_op = 0*vd_op;
id_op = (2^0.5)*real(    (Ptt + 1j*Qtt)./(Vmagtt.*exp(1j*phtt) ) ); % Id at operating point
iq_op = -(2^0.5)*imag( (Ptt + 1j*Qtt)./(Vmagtt.*exp(1j*phtt) ) ); % Iq at operating point
i0_op = 0*id_op;
% Approximate the large-signal step response by adding
% the small-signal values to the operating point.
vds = vd + vd_op*ones(1,length(t));
vqs = vq + vq_op*ones(1,length(t));
v0s = v0 + v0_op*ones(1,length(t));
ids = id + id_op*ones(1,length(t));
iqs = iq + iq_op*ones(1,length(t));
i0s = i0 + i0_op*ones(1,length(t));
ps = 0.5*(vds.*ids + vqs.*iqs + 2*v0s.*i0s); % active powers of generators


%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Visualize results:
%%%%%%%%%%%%%%%%%%%%%%%%
% Plot step response (active power of generators)
set(0,'defaulttextinterpreter','latex')
set(0,'defaultfigurecolor',[1 1 1])
set(0,'defaultaxesfontsize',9);
figure(2);
plot(t,ps/baseVA,'LineWidth',0.5);
xlabel('Time [s]');
ylabel('Active powers of generators [pu]');

axesHandles = findall(0,'type','axes');
set(axesHandles,'TickLabelInterpreter', 'latex')
