function mpc_out = line2shunt(mpc_in)
% LINE2SHUNT operates on a MatPower database.
% It converts capacitors that appear on branches to shunt elements
% on the adjacent buses.

% Define standard column numbers used in the MatPower database
F_BUS = 1; T_BUS = 2; BR_B = 5; BUS_I = 1; BS = 6;

mpc_out = mpc_in;
for kk=1:size(mpc_in.branch,1)
    fb = mpc_in.branch(kk,F_BUS); % from bus
    tb = mpc_in.branch(kk,T_BUS); % to bus
    bsp = mpc_in.branch(kk,BR_B)*(mpc_in.baseMVA)/2; % additional shunt capacitance in MVAr injected at 1p.u.
    mpc_out.branch(kk,BR_B) = 0; % eliminate branch capacitance
% Find buses on bus list:
    bus_ind1 = find( mpc_in.bus(:,BUS_I) == fb );
    bus_ind2 = find( mpc_in.bus(:,BUS_I) == tb );  
% Add shunt capacitors to buses:
    mpc_out.bus(bus_ind1,BS) = mpc_out.bus(bus_ind1,BS) + bsp;
    mpc_out.bus(bus_ind2,BS) = mpc_out.bus(bus_ind2,BS) + bsp;
end

end
