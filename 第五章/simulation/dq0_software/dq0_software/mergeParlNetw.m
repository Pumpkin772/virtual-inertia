function [A, B, C, D] = mergeParlNetw(A1, B1, C1, D1, A2, B2, C2, D2)
% MERGEPARLNETW Merges two networks connected in parallel.
% 
% Usage:
%       [A, B, C, D] = MERGEPARLNETW(A1, B1, C1, D1, A2, B2, C2, D2)
% 
% where
%       A1, B1, C1, D1 - system matrices of network 1
%       A2, B2, C2, D2 - system matrices of network 2
%
% Outputs:
%       A, B, C, D - system matrices of merged network
% 
% 
% **** Comments: ****
% The networks are described by state-space models.
% The model inputs (voltages in dq0 coordinates) are the same for the two networks.
% The model outputs (injected currents in dq0 coordinates) are summed entry-wise.
% Merged networks cannot be used in functions 'createQS', and 'createYbus'


narginchk(8,8)

A = [A1 sparse(size(A1,1),size(A2,2)); ...
    sparse(size(A2,1),size(A1,2)) A2];
B = [B1; B2];
C = [C1 C2];
D = D1+D2;

end
