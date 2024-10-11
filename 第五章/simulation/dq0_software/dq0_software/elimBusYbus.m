function YbusR = elimBusYbus(Ybus, subset)
% ELIMBUSYBUS Eliminates disconnected buses in the system admittance matrix Ybus.
%
% Usage:
%       YbusR = elimBusYbus(Ybus, subset)
% 
% where
%       Ybus - original admittance matrix
%       subset - a vector of bus indices specifying the buses
%           to be included in the new model. Voltages and currents
%           of buses in this list are NOT eliminated, and will be
%           the inputs and outputs of the new model
%
% Outputs:
%       YbusR - the new admittance matrix
% 
% 
% **** Comments: ****
% The original admittance matrix describes a quasi-static model and
% relates bus voltages to bus currents by
% I = Ybus*V
% where I and V are vectors of length N (number of buses).
% The elements of V and I are the bus voltage phasors
% and the bus current phasors, respectively.
%
% After bus elimination the new relation is
% I = YbusR*V
% where I and V are new vectors of length M, and
% M is the number of remaining buses (equal to the length of 'subset')


narginchk(2,2)

N = size(Ybus,1);
tmp = ones(N,1);
tmp(subset) = 0;
elt_ind = find(tmp);
UU = speye(N);
UU(:,elt_ind)=[];
QQ = speye(N);
QQ(:,subset)=[];
YbusR = (UU')*Ybus*UU - ((UU')*Ybus*QQ)*...
    inv((QQ')*Ybus*QQ)*((QQ')*Ybus*UU);

end
