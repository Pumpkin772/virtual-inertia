function Ybus = createYbus(A, B, C, D)
% CREATEYBUS Creates admittance matrix Ybus based on system matrices of
% the full dq0 model.
% 
% Usage:
%       Ybus = createYbus(A, B, C, D)
% 
% where
%       A, B, C, D - system matrices of dq0 model, computed by functions
%           'ssNetw' or 'ssNetwSym'. These matrices can be full, sparse,
%           numeric, or symbolic
%
% Outputs:
%       Ybus - the network Ybus matrix that is numeric or symbolic,
%           depending on the inputs
% 
% See also: ssNetw, ssNetwSym


narginchk(4,4)

N = size(B,2)/3; % number of buses
MM = size(A,1); % number of states

% Verify legal structure of the input matrices
valid = 1;
valid = valid & isempty(find(A(1:(MM*2/3),(MM*2/3+1):end), 1));
valid = valid & isempty(find(A((MM*2/3+1):end,1:(MM*2/3)), 1));
valid = valid & isempty(find(B((MM*2/3+1):end,1:(2*N)), 1));
valid = valid & isempty(find(C(1:(2*N),(MM*2/3+1):end), 1));
if (~valid)
    error('createYbus:InvalidStructure', ...
        'Invalid structure of input matrices.');
end

% Remove parts associated with the zero component
Aqs = A(1:(MM*2/3),1:(MM*2/3));
Bqs = B(1:(MM*2/3),1:(2*N));
Cqs = C(1:(2*N),1:(MM*2/3));
Dqs = D(1:(2*N),1:(2*N));

P1 = -Cqs*(Aqs\Bqs)+Dqs;
if (isa(P1,'sym'))
    P1 = simplify(P1);
end
% Verify P1 structure
if (isa(P1,'sym'))
    u1 = isequal( P1(1:N,1:N), P1((N+1):2*N,(N+1):2*N) );
    u2 = isequal( P1((N+1):2*N,1:N) , -P1(1:N,(N+1):2*N));
else
    eps = 1e-7;
    e1 = max(max(abs(P1(1:N,1:N) - P1((N+1):2*N,(N+1):2*N))));
    e2 = max(max(abs(P1((N+1):2*N,1:N) + P1(1:N,(N+1):2*N))));
    e1a = max(max(abs( P1(1:N,1:N))));
    e1b = max(max(abs( P1((N+1):2*N,(N+1):2*N))));
    e2a = max(max(abs( P1((N+1):2*N,1:N))));
    e2b = max(max(abs( P1(1:N,(N+1):2*N))));
    if (max(e1a,e1b) == 0)
        u1 = 1;
    else
        u1 = (e1/max(e1a,e1b)) < eps;
    end
    if (max(e2a,e2b) == 0)
        u2 = 1;
    else
        u2 = (e2/max(e2a,e2b)) < eps;
    end
end
if (~(u1&&u2))
    error('createYbus:InvalidStructure', ...
        'Invalid structure of matrix P1.');
end
Ybus = P1(1:N,1:N) + 1j*P1((N+1):2*N,1:N);
if (isa(Ybus,'sym'))
    Ybus = simplify(Ybus);
end

end
