function mpcOut = switchBusInd(mpcIn, bus1, bus2)
% SWITCHBUSIND Switches between two buses in a MatPower database.
% 
% Usage:
%       mpcOut = SWITCHBUSIND(mpcIn, bus1, bus2)
% where
%       mpcIn - MatPower database
%       bus1, bus2 - numbers of buses to switch
% 
% Outputs:
%       mpcOut - new MatPower database
% 
% 
% **** Comments: ****
% Useful for numbering the reference bus (slack bus) as 1

narginchk(3,3)

% Define standard column numbers used in the MatPower database
BUS_I = 1; GEN_BUS = 1; F_BUS = 1; T_BUS = 2; BR_B = 5; BUS_I = 1; BS = 6;

mpcOut = mpcIn;

ind1 = find(mpcIn.bus(:,BUS_I) == bus1);
ind2 = find(mpcIn.bus(:,BUS_I) == bus2);
mpcOut.bus(ind1,BUS_I) = bus2;
mpcOut.bus(ind2,BUS_I) = bus1;

ind1 = find(mpcIn.gen(:,GEN_BUS) == bus1);
ind2 = find(mpcIn.gen(:,GEN_BUS) == bus2);
mpcOut.gen(ind1,GEN_BUS) = bus2;
mpcOut.gen(ind2,GEN_BUS) = bus1;

ind1 = find(mpcIn.branch(:,F_BUS) == bus1);
ind2 = find(mpcIn.branch(:,F_BUS) == bus2);
mpcOut.branch(ind1,F_BUS) = bus2;
mpcOut.branch(ind2,F_BUS) = bus1;

ind1 = find(mpcIn.branch(:,T_BUS) == bus1);
ind2 = find(mpcIn.branch(:,T_BUS) == bus2);
mpcOut.branch(ind1,T_BUS) = bus2;
mpcOut.branch(ind2,T_BUS) = bus1;

if isfield(mpc,'dcline')
    ind1 = find(mpcIn.dcline(:,F_BUS) == bus1);
    ind2 = find(mpcIn.dcline(:,F_BUS) == bus2);
    mpcOut.dcline(ind1,F_BUS)=bus2;
    mpcOut.dcline(ind2,F_BUS)=bus1;

    ind1 = find(mpcIn.dcline(:,T_BUS) == bus1);
    ind2 = find(mpcIn.dcline(:,T_BUS) == bus2);
    mpcOut.dcline(ind1,T_BUS)=bus2;
    mpcOut.dcline(ind2,T_BUS)=bus1;
end

end
