% This script demonstrates the use of several basic functions:
% - ssNetw
% - ssNetwSym
% - createYbus
% - createQS
% - elimBus
% - elimBusesYbus
% - longLine

clc;

% Choose an example:
choose_example = 1;

if (choose_example==1)
% Example: create a dq0 state-space model of a 3 bus system.
% branch 1->2 contains inductance, resistance, and a transformer.
% branch 1->3 does not exist
% branch 2->3 does not exist
% bus 1 has a shunt resistive-capacitive element
% bus 2 has a shunt inductance
% bus 3 has a shunt resistance

% Define network parameters :
% (See documentation in 'ssNetw')  
    ws = 2*pi*60; % network frequency at operating point

% Branch elements
    R12 = 3; % resistance of branch 1->2
    R23 = 0;
    R13 = 0;
    L12 = 0.1; % inductance of branch 1->2
    L13 = inf;  % branch does not exist    
    L23 = 0;  % branch does not exist (0 is like inf)
    Rb = [0 R12 R13; 0 0 R23; 0 0 0]; % only upper diagonal is defined
    Lb = [0 L12 L13; 0 0 L23; 0 0 0];
    Tau = zeros(3);
    Tau(2,1) = 0.9; % a transformer at branch 1->2 on bus 2
    Tau(1,3) = 1; % unity transfer ratios are ignored

% Shunt elements
% bus 1 has a shunt resistive-capacitive element
% bus 2 has a shunt inductance
% bus 3 has a shunt resistance
    R_bus = [0 0 6];
    L_bus = [0 5 0];
    C_bus = [2 inf 0];  % 0 and inf are the same
    Rtil_bus = [7 0 0];
    
% Build minimal state-space model
    [A, B, C, D] = ssNetw(R_bus, L_bus, C_bus, Rtil_bus, Rb, Lb, Tau, ws);
         
% Now eliminate the 1st and 2nd buses.
% This produces an equivalent model in which
% the currents id(1), iq(1), i0(1) and id(2), iq(2), i0(2) are zeroed.
% The inputs and outputs of the new model correspond to the 3rd bus.
% (for more details see documentation in 'elimBus')
    subset = [3];  % remaining buses
    [Ar, Br, Cr, Dr ] = elimBus(A, B, C, D, subset);
    sys = ss(full(Ar), full(Br), full(Cr), full(Dr)); % create state-space model
    sys = reduce(sys,'MaxError',0); % reduce model, eliminate redundant states
    sys
    
% The result: a static gain model.
% This model describes a resistor to ground, which is the shunt element of bus 3.
end

if (choose_example==2) 
% Example: Symbolic state-space model.
% Single bus with shunt inductance and resistance.
% No branches. quasi-static model (Ybus matrix)
% is computed from the state-space matrices.

    syms R L
    
    R_bus = [R]; % shunt resistance
    L_bus = [L]; % shunt inductance
    C_bus = [0]; % infinite shunt capacitance 
    Rtil_bus = [0]; % infinite shunt series resistance
    Rb = 0;
    Lb = 0;
    Tau = 0;

% Generate state-space model:
    [Asym, Bsym, Csym, Dsym] = ssNetwSym(R_bus, L_bus, C_bus, Rtil_bus, Rb, Lb, Tau);
    Asym
    Bsym
    Csym
    Dsym
      
% Create quasi-static model
    [Aqs, Bqs, Cqs, Dqs] = createQS(Asym, Bsym, Csym, Dsym);
    quasi_static_gain = Dqs
    
% Extract Ybus matrix from state-space model:
    Ybus = createYbus(Asym, Bsym, Csym, Dsym);
    Ybus
end

if (choose_example==3)
% Example:  Symbolic State-Space Model
% Two bus network.
% No shunt elements.
% branch 1->2 contain an inductor and a transformer.
    
    syms L tr
   
% No shunt elements
    R_bus = [0 0];
    L_bus = [0 0];
    C_bus = [0 0];
    Rtil_bus = [0 0];
    
% Branches
    Rb = [0 0; 0 0];
    Lb = [0 L; 0 0];
    Tau = [0 tr; 0 0]; % transformer on branch 1->2, located on bus 1
       
% Generate state-space model
    [Asym, Bsym, Csym, Dsym] = ssNetwSym(R_bus, L_bus, C_bus, Rtil_bus, Rb, Lb, Tau);
   
% Display results
    disp('Symbolic state-space model:')
    Asym
    Bsym
    Csym
    Dsym
    
% Extract Ybus matrix from dq0 state-space model:
    Ybus = createYbus(Asym,Bsym,Csym,Dsym);
    Ybus
% Evaluate Ybus with no transformer (tr=1)
% (producing the well-known Ybus structure)
    Ybus_no_transformer = subs(Ybus,tr,1);
    Ybus_no_transformer
end

if (choose_example==4)
% Example: Symbolic model with bus elimination.
% The network includes 3 buses, and the 2nd bus is eliminated.
% The objective is to find an equivalent state space model
% in which the currents in bus 2 are zeroed, so the
% reduced model includes only two buses (the original bus 1, bus 3).

% The (original) network consists of 3 buses:
% buses 1 & 2 are connected by an inductance L
% buses 2 & 3 are connected by an inductance L
% bus 2 include a shunt resistance to ground R2

    syms R L
  
% Define network parameters
    R_bus = [0 R 0].';
    L_bus = [0 0 0];
    C_bus = [0 0 0];
    Rtil_bus = [0 0 0];
    Rb = [0 0 0; 0 0 0; 0 0 0];
    Lb = [0 L 0; 0 0 L; 0 0 0];
    Tau = zeros(3);
      
% Build minimal state-space model (original one)
    [Asym, Bsym, Csym, Dsym] = ssNetwSym(R_bus, L_bus, C_bus, Rtil_bus, Rb, Lb, Tau);
    
% Eliminate the 2nd bus
% This produces an equivalent model in which
% the currents id(2), iq(2), i0(2) are zeroed.
% The inputs and outputs of the new model
% correspond to the 1st and 3rd buses.
% (for more details see documentation in 'elimBus')
    subset = [1 3];
    [Asym_eq, Bsym_eq, Csym_eq, Dsym_eq] = elimBus(Asym, Bsym, Csym, Dsym, subset);
         
% Additional transformation to balance '2' and '1/2'
% entries in matrices B and C.
% (This is just cosmetics and not really significant)
    TT = sym([1/2 0 0 0 0 0; 
              0 1/2 0 0 0 0;
              0 0 1/2 0 0 0;
              0 0 0 1/2 0 0;
              0 0 0 0   1 0;
              0 0 0 0   0 1;]);
    Asym_eq = inv(TT)*Asym_eq*TT;
    Bsym_eq = inv(TT)*Bsym_eq;
    Csym_eq = Csym_eq*TT;
    
    Asym_eq
    Bsym_eq
    Csym_eq
    Dsym_eq   
      
% Extract Ybus matrix from dq0 state-space model:
    Ybus = createYbus(Asym_eq, Bsym_eq, Csym_eq, Dsym_eq);

% The analytic Ybus in this example is given as follows.
% Note that this matrix has a non-traditional
% structure since one bus has been reduced.
    syms ws
    Ybus_analytic = [1j*ws*L+R -R; -R 1j*ws*L+R]/(2*1j*ws*L*R-ws^2*L^2);
    
% Test equality of Ybus matrices by substitution of random parameters:
    L = sym(1/round(1+100*rand(1,1)));
    R = sym(round(1+100*rand(1,1)));
    ws = 2*pi*55;
% Evaluate the difference in Ybus matrices:
    error_in_ybus = max(max(abs(double(subs(Ybus-Ybus_analytic)))));
    error_in_ybus
    disp('Should be very small.');
end

if (choose_example==5)
% Example: Symbolic model with bus elimination.
% In this example one bus is eliminated from the model,
% and the corresponding admittance matrix Ybus is
% computed in two different methods:
% a) by eliminating buses in the dynamic model
% and then computing Ybus of this new model.
% b) by computing Ybus of the original dynamic model
% and eliminating it directly.
% The results are shown to be equal
    
    syms L C Rtil
  
% Define inputs
    R_bus = [0 0];
    L_bus = [0 0];
    C_bus = [0 C];
    Rtil_bus = [0 Rtil];
    Rb = [0 0; 0 0];
    Lb = [0 L; 0 0];
    Tau = zeros(2);
                  
% Build state-space model (original one, before bus elimination)
    [Asym, Bsym, Csym, Dsym] = ssNetwSym(R_bus, L_bus, C_bus, Rtil_bus, Rb, Lb, Tau);
              
% Eliminate the 2nd bus in the dynamic model
    subset = 1; % bus 1 remains, bus 2 is eliminated
    [Ar, Br, Cr, Dr] = elimBus(Asym, Bsym, Csym, Dsym, subset);
      
% Create Ybus of eliminated model
    Ybus_el = createYbus(Ar, Br, Cr, Dr);
    
% Create original Ybus and eliminate it directly
    Ybus_ppp = createYbus(Asym, Bsym, Csym, Dsym);
    Ybus = elimBusYbus(Ybus_ppp, subset);
    
% Compare the results
% (Although it's hard to see, these expressions are equal)
    Ybus
    Ybus_el
       
% Test equivalence of the result by substituting some random numbers:  
    Rtil = 1+100*rand(1);
    L = 1+100*rand(1);
    C = 1+100*rand(1);
    ws = 2*pi*50;
      
    result1 = double(subs(Ybus));
    result2 = double(subs(Ybus_el));
    error = abs(result1 - result2);
    error  % should be zero
end

if (choose_example==6) 
% Example: Testing speed and memory with large networks.
% All values are random.
% N branches are connected in random (the others are not connected)
% For 3000 bus system (N=3000):
% Estimated computing time is ~4-8 min.
% Max memory usage during computation is about 6 GB
    
    N = 500;  % number of buses
    
    R_bus = 1+rand(N,1);
    C_bus = 1+rand(N,1);
    L_bus = 1+rand(N,1);
    Rtil_bus = 1+rand(N,1);
    Rb = sparse(N,N);
    Lb = sparse(N,N);
    Tau = sparse(N,N);
    ws = [2*pi*50];
% Connect N branches randomly:
    for ii=1:N
        ind1 = ceil(N*rand(1));
        ind2 = ceil(N*rand(1));
        if (ind2<ind1)  % put elements in the upper diagonal
            Lb(ind2,ind1) = 1+rand(1,1);
            Rb(ind2,ind1) = 1+rand(1,1);          
        elseif (ind1<ind2)
            Lb(ind1,ind2) = 1+rand(1,1);
            Rb(ind1,ind2) = 1+rand(1,1);
        end
        if (ind1 ~= ind2)
            Tau(ind1,ind2) = 0.5+rand(1,1);
        end
    end
    
    tic
% Generate state-space model:
    [A, B, C, D] =  ssNetw(R_bus, L_bus, C_bus, Rtil_bus, Rb, Lb, Tau, ws);
% Extract Ybus matrix from state-space model:
    Ybus = createYbus(A, B, C, D);
    disp('done');
    toc
      
% Plot sparsity patterns:
    figure(1);
    subplot(2,2,1); spy(A); title('A');
    subplot(2,2,2); spy(B); title('B');
    subplot(2,2,3); spy(C); title('C');
    subplot(2,2,4); spy(D); title('D');
end

if (choose_example==7) 
% Additional testing for function 'elimBus'.
% testing the condition where eliminated buses have
% zero values of the diagonal of D. In this case
% elimination is achieved through LU decomposition
% of the constraints (see function 'eliminate_buses')

% Define a network with N buses, N-1 braches, and no loads.
% Each pair of buses [n-1,n] are connected with an inductor L
% (The network is simply several inductors connected in series)

% The objective in this example is to eliminate the middle buses.
% The resulting network has two buses that are connected through an inductor
% (N-1)*L
    
    N = 101; % number of buses
   
% Define network
    R_bus = zeros(N,1);
    L_bus = zeros(N,1);
    C_bus = zeros(N,1);
    Rtil_bus = zeros(N,1);
    Rb = sparse(N,N);
    Lb = sparse(N,N);
    L = 1;
    for n=2:N
        Lb(n-1,n) = L;
    end
    Tau = sparse(N,N);
    ws = 2*pi*50;
    [A1, B1, C1, D1] = ssNetw(R_bus, L_bus, C_bus, Rtil_bus, Rb, Lb, Tau, ws);
    
% Now eliminate buses 2...N-1
    subset = [1,N]; % buses to keep
    if (N==1)
        subset = 1;
    end
    [A2, B2, C2, D2] = elimBus(A1, B1, C1, D1, subset);
    
% Define 2 bus network in which buses 1,2 are connected by
% an indcutor of (N-1)*L
    R_bus = [0 0];
    L_bus = [0 0];
    C_bus = [0 0];
    Rtil_bus = [0 0];
    Rb = sparse(2,2);
    Lb = sparse(2,2);
    Lb(1,2)=(N-1)*L;
    Tau = sparse(2,2);
    ws = 2*pi*50;
    [A3, B3, C3, D3] = ssNetw(R_bus, L_bus, C_bus, Rtil_bus, Rb, Lb, Tau, ws);
       
% State space model A2, B2, C2, D2 is identical to A3, B3, C3, D3
    clc
    full(A3)
    full(A2)
    full(B3)
    full(B2)
    full(C3)
    full(C2)
    
    size_of_original_A_matrix = size(A1)
    size_of_new_A_matrix = size(A2)
end

if (choose_example==8) 
% This is a dq0 model of a long transmission line.
% The line is constructed from Pi sections
% with shunt capacitors C/2 on both sides,
% and inductance L connected in series.

% The line is constructed as a network
% with N buses (N-1 Pi sections).

% The middle buses are eliminated to obtain
% an equivalent model of the long transmission (with two buses).
               
    Npi = 50; % number of Pi sections
    N = Npi+1; % number of buses
    
    L = 1e-4/Npi; % [H] inductance of Pi section
    C = 1e-4/Npi; % [F] total capacitance of Pi section (each side C/2)
          
% Expected time delay with infinite number of Pi sections:
    expected_delay = Npi*(L*C)^0.5; 
% Characteristic impedance:
    Zc = (L/C)^0.5; % [ohm]
% Series resistance of shunt capacitors
    Rser = Zc/20;  % small value
      
% Define network
% Transmission line modeled by Pi sections
    R_bus = zeros(N,1);
    L_bus = zeros(N,1);
    C_bus = C*ones(N,1);
    C_bus(1) = C/2; C_bus(end) = C/2;
    Rtil_bus = 2*Rser*ones(N,1);
    Rtil_bus(1) = Rser; Rtil_bus(end) = Rser;
    Rb = sparse(N,N);
    Lb = sparse(N,N);
    for n=2:N
        Lb(n-1,n)=L;
    end
    Tau = sparse(N,N);
    ws = 2*pi*50;
    [A1, B1, C1, D1] = ssNetw(R_bus, L_bus, C_bus, Rtil_bus, Rb, Lb, Tau, ws);
    
% Now eliminate buses 2...N-1
    subset = [1,N]; % buses to keep
    if (N==1)
        subset = 1;
    end
    [A2, B2, C2, D2] = elimBus(A1, B1, C1, D1, subset);
    
% The network is terminated with a resistor
% Zc = (L/C)^0.5 (the characteristic impedance)
% This resistor creates a feedback from the
% currents of the 2nd bus to the voltage of the 2nd bus.
% Here is the resulting dynamic system including this feedback:
% d/dt x = A2*x + B2*u
% y = C2*x + D2*u
% Feedback:  u = W1*Vin - W2*Zc*W3*y,
% where: Vin = [v_{d,1}; v_{q,1}; v_{0,1}],
% Zc characteristic impedance
% W1, W2, W3 matrices that selects proper elements in the vectors.
% Algebra:
% y = C2*x+D2*(W1*Vin - W2*Zc*W3*y)
% y = C2*x+D2*W1*Vin - D2*W2*Zc*W3*y
% (eye+D2*W2*Zc*W3)*y = C2*x + D2*W1*Vin
% define P = eye+D2*W2*Zc*W3:
% P*y = C2*x + D2*W1*Vin
% y = invP*C2*x + invP*D2*W1*Vin
% u = -W2*Zc*W3*invP*C2*x + (W1-W2*Zc*W3*invP*D2*W1)*Vin
% dx/dt = (A2-B2*W2*Zc*W3*invP*C2)*x + B2*(W1-W2*Zc*W3*invP*D2*W1)*Vin
% define W4 signal routing matrix that selects the output
% to be dq0 currents of the 2nd bus:
% Iout = [i_{d,2}; i_{q,2}; i_{0,2}]
% Iout = W4*y
% Iout = W4*invP*C2*x + W4*invP*D2*W1*Vin
% New state-space model with feedback
% d/dt x = Af*x + Bf*Vin
% Iout = Cf*x + Df*Vin    
% Af = A2-B2*W2*Zc*W3*invP*C2
% Bf = B2*(W1-W2*Zc*W3*invP*D2*W1)
% Cf = W4*invP*C2
% Df = W4*invP*D2*W1 

% Define signal routing matrices
    W1 = zeros(6,3);
    W1(1,1)=1; W1(3,2)=1; W1(5,3)=1;
    W2 = zeros(6,3);
    W2(2,1)=1; W2(4,2)=1; W2(6,3)=1;
    W3=W2';
    W4=-W3;
       
    P = eye(size(D2))+D2*W2*Zc*W3;
    invP = inv(P);
    
% New state-space with feedback
% input vector Vin = [v_{d,1}; v_{q,1}; v_{0,1}]
% output vector Iout = [i_{d,2}; i_{q,2}; i_{0,2}]
    Af = A2-B2*W2*Zc*W3*invP*C2;
    Bf = B2*(W1-W2*Zc*W3*invP*D2*W1);
    Cf = W4*invP*C2;
    Df = W4*invP*D2*W1;
    
    sys = ss(Af,Bf,Cf,Df);
    sys = reduce(sys,'MaxError',0.1);
    sys_d = sys(1,1); % from v_{d,1} to i_{d,2}
    sys_q = sys(2,2); % from v_{q,1} to i_{q,2}
    sys_0 = sys(3,3); % from v_{0,1} to i_{0,2}
    sys_qd = sys(2,1); % from v_{d,1} to i_{q,2}
    
    set(0,'defaulttextinterpreter','latex')
    set(0,'defaultfigurecolor',[1 1 1])
    set(0,'defaultaxesfontsize',9);
    figure(1);
    step(sys_d,expected_delay*5);
    set(gca,'TickLabelInterpreter', 'latex')
    
    number_of_pi_section = Npi
    expected_delay
end

if (choose_example==9) 
% Implements a dynamic model of a long transmission line
    Rx = 2e-4; % [Ohm/m]
    Lx = 1e-6; % [H/m]
    Cx = 1e-11; % [F/m]
    len = 5e5; % [m]
    N = 10;
    ws = 2*pi*50;
    
    R = Rx*len/N;
    L = Lx*len/N;
    C = Cx*len/N;
    
    expected_delay = len*(Lx*Cx)^0.5; % characteristic impedance: 
    Rc = (Lx/Cx)^0.5; % [Ohm]
    Rc = Rc;  % this value can be changed to obtain travelling waves
    
% Transmission line state-space model
    [At, Bt, Ct, Dt] = longLine(Rx, Lx, Cx, len, N, ws);
    
% With feedback of load resistor Rc:
% input is [v_{d,S}; v_{q,S}; v_{0,S}];
% output is [i_{d,R}; i_{q,R}; i_{0,R}];
    ey = speye(3);
    zr = sparse(3,3);
    Af = At - Rc*Bt*[zr , zr ; zr, ey]*Ct;
    Bf = Bt*[ey ; zr];
    Cf = [zr , ey]*Ct;
    Df = sparse(3,3);
    
% Compute 1kV step response for step in v_{d,S}
    sys = ss(full(Af), full(Bf), full(Cf), full(Df));
    [ir, t] = step(sys(:,1),0.02);
    
    set(0,'defaulttextinterpreter','latex')
    set(0,'defaultfigurecolor',[1 1 1])
    set(0,'defaultaxesfontsize',9);
    figure(1);
    plot(t,1000*ir,'-k','Color',[0 0 0],'LineWidth',0.5);
    xlabel('Time [s]')
    ylabel('Receiving-end $dq0$ current [A]')
    set(gca,'TickLabelInterpreter', 'latex')
end
