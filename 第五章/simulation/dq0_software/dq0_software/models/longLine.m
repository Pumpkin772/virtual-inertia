function [A, B, C, D] = longLine(Rx, Lx, Cx, len, N, ws)
% LONGLINE Implements an approximate (lumped) dynamic model of a long power line.
%
% Usage:
%       [A, B, C, D] = LONGLINE(Rx, Lx, Cx, len, N, ws)
% 
% where
%       Rx - resistance per unit length [Ohm/m]
%       Lx - inductance per unit length [H/m]
%       Cx - capacitance per unit length [F/m]
%       len - transmission line length [m]
%       N - number of T sections 
%       ws - nominal system frequency [rad/s]
%
% Outputs:
%       A, B, C, D - system matrices
% 
% See also: ssNetw


narginchk(6,6)

R = Rx*len/N;
L = Lx*len/N;
C = Cx*len/N;
Ap = sparse([-R/L ws 0; -ws -R/L 0; 0 0 -R/L]);
App = sparse([0 ws 0; -ws 0 0; 0 0 0]);
I3 = speye(3);

A = sparse(6*N+3,6*N+3);
A(1:3,1:3) = Ap;
A(1:3,4:6) = -(2/L)*I3;
A((6*N+1):(6*N+3),(6*N-2):(6*N)) = (2/L)*I3;
A((6*N+1):(6*N+3),(6*N+1):(6*N+3)) = Ap;
for kk=1:N,
    A((6*kk-2):(6*kk),(6*kk-5):(6*kk-3)) = (1/C)*I3;
    A((6*kk-2):(6*kk),(6*kk-2):(6*kk)) = App;
    A((6*kk-2):(6*kk),(6*kk+1):(6*kk+3)) = -(1/C)*I3;
end
for kk=2:N,
    A((6*kk-5):(6*kk-3),(6*kk-8):(6*kk-6)) = (1/L)*I3;
    A((6*kk-5):(6*kk-3),(6*kk-5):(6*kk-3)) = Ap;
    A((6*kk-5):(6*kk-3),(6*kk-2):(6*kk)) = -(1/L)*I3;
end

B = sparse(6*N+3,6);
B(1:3,1:3) =  (2/L)*I3;
B((6*N+1):(6*N+3),4:6) = -(2/L)*I3;
C = sparse(6,6*N+3);
C(1:3,1:3) = I3;
C(4:6,(6*N+1):(6*N+3)) = -I3;
D = sparse(6,6);

end
