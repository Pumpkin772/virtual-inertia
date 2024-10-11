function [Aphi, Bphi, Cphi, Dphi, Gphi] = infBus()
% INFBUS Constructs a linear model of an infinite bus (constant voltage source).
%
% Usage:
%       [A, B, C, D, G] = infBus()
% 
% Outputs:
%       A, B, C, D, G - system matrices


Aphi = sparse(0,0);
Bphi = sparse(0,3);
Cphi = sparse(3,0);
Dphi = sparse(3,3);
Gphi = sparse(0,0);

end
