% Examples based on the MatPower database.
% This file loads Matpower test-case networks
% and constructs their dynamic models, in the dq0 reference frame.
%
% This script requires the MatPower package to be installed (available online).

clear all;
define_constants; % MatPower function. Defines MatPower constants.

% Choose MatPower test-case network
switch_case = 8;

switch switch_case
    case 1
        filename = 'case4gs';
        mpc = loadcase(filename);
    case 2
        filename = 'case5';
        mpc = loadcase(filename);
    case 3
        filename = 'case6ww';
        mpc = loadcase(filename);
    case 4
        filename = 'case9';
        mpc = loadcase(filename);
    case 5
        filename = 'case9Q';
        mpc = loadcase(filename);
    case 6
        filename = 'case9target';
        mpc = loadcase(filename);
    case 7
        filename = 'case14';
        mpc = loadcase(filename);
    case 8
        filename = 'case24_ieee_rts';
        mpc = loadcase(filename);
        % remove duplicate branches
        mpc.branch(26,:)=[];
        mpc.branch(32,:)=[];
        mpc.branch(33,:)=[];
        mpc.branch(34,:)=[];
    case 9
        filename = 'case_ieee30';
        mpc = loadcase(filename);
    case 10
        filename = 'case30Q';
        mpc = loadcase(filename);
    case 11
        filename = 'case39';
        mpc = loadcase(filename);
    case 12
        filename = 'case57';
        mpc = loadcase(filename);
        % remove duplicate branches
        mpc.branch(20,:)=[];
        mpc.branch(35,:)=[];
    case 13
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
    case 14
        filename = 'case2383wp';
        % This example may take several minutes to compute (~2 min.)
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
    case 15
        filename = 'case2736sp';
        % This example may take several minutes to compute (~3.5 min.),
        % and requires at least 4 GB of free memory.
        mpc = loadcase(filename);
        % remove duplicate branches
        mpc.branch(209,:)=[];
        mpc.branch(233,:)=[];
        mpc.branch(920,:)=[];
        mpc.branch(1035,:)=[];
        mpc.branch(2800,:)=[];
        mpc.branch(1394,:)=[];
        mpc.branch(2542,:)=[];
        mpc.branch(2662,:)=[];
        mpc.branch(2630,:)=[];
        % remove phase shifts in transformers:
        mpc.branch(:,SHIFT)=0;        
    otherwise
        disp('Unknown Matpower test-case');
        return
end

fs = 50;  % [Hz] nominal system frequency
ws = 2*pi*fs; % [rad/s]

% Compute dq0 state-space model
[A, B, C, D, YbusPU] = ssNetwMatPower(mpc, ws);
% [Aqs, Bqs, Cqs, Dqs] = createQS(A, B, C, D);
% nnz(Dqs)/(size(Dqs,1)*size(Dqs,2))*100
% nnz(Dqs)


% % Compute admittance matrix (Ybus) using Matpower function
% [Ybus, Yf, Yt] = makeYbus(mpc);
% 
% % Compare results for admittance matrices:
% % The admittance matrix YbusPU (generated above) should be identical
% % to the admittance matrix generated directly by MatPower.
% % If the state-space model is correct then these matrices should be nearly equal.
% Ybus_diff = (Ybus-YbusPU) / max(max(abs(Ybus))) ;
% Ybus_error = max(abs(Ybus_diff(:)));
% 
% filename
% full(Ybus_error)
% % disp('if the dynamic model is correct');
% % disp('the normalized error should be small.');
% % disp('(typically smaller than 1e-6)');
