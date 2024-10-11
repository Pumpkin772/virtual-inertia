function [Aqs, Bqs, Cqs, Dqs] = createQS(A, B, C, D)
% CREATEQS creates a quasi-static model from a dynamic dq0 model.
%
% Usage:
%       [Aqs, Bqs, Cqs, Dqs] = create_quasi_static(A, B, C, D)
% 
% where
%       A, B, C, D - system matrices of dq0 model, computed by functions
%           'ssNetw' or 'ssNetwSym'. These matrices can be full, sparse,
%           numeric, or symbolic
% 
% Outputs:
%       Aqs, Bqs, Cqs, Dqs - system matrices of the quasi-static model that
%       are numeric or symbolic, depending on the inputs
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
    error('createQS:InvalidStructure', ...
        'Invalid structure of input matrices.');
end

% Remove parts associated with the zero component
At = A(1:(MM*2/3),1:(MM*2/3));
Bt = B(1:(MM*2/3),1:(2*N));
Ct = C(1:(2*N),1:(MM*2/3));
Dt = D(1:(2*N),1:(2*N));

Aqs = [];
Bqs = [];
Cqs = [];
Dqs = -Ct*(At\Bt)+Dt;
if (isa(Dqs,'sym'))
    Dqs = simplify(Dqs);
end
% Verify Dqs structure
if (isa(Dqs,'sym'))
    u1 = isequal(Dqs(1:N,1:N), Dqs((N+1):2*N,(N+1):2*N));
    u2 = isequal(Dqs((N+1):2*N,1:N), -Dqs(1:N,(N+1):2*N));
else
    eps = 1e-9;
    e1 = max(max(abs(Dqs(1:N,1:N) - Dqs((N+1):2*N,(N+1):2*N))));
    e2 = max(max(abs(Dqs((N+1):2*N,1:N) + Dqs(1:N,(N+1):2*N))));
    e1a = max(max(abs(Dqs(1:N,1:N))));
    e1b = max(max(abs(Dqs((N+1):2*N,(N+1):2*N))));
    e2a = max(max(abs(Dqs((N+1):2*N,1:N))));
    e2b = max(max(abs(Dqs(1:N,(N+1):2*N))));
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
    error('createQS:InvalidStructure', ...
        'Invalid structure for matrix Dqs.');   
end

end
