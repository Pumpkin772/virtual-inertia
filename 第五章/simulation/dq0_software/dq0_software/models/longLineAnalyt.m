function sys = longLineAnalyt(Rx, Lx, Cx, len, ws, w)
% LONGLINEANALYT Creates frequency-response data model of
% a long transmission line in the dq0 reference frame.
%
% Usage:
%       sys = longLineAnalyt(Rx, Cx, Lx, l, ws, w)
% 
% where
%       Rx - resistance per unit length [Ohm/m]
%       Lx - inductance per unit length [H/m]
%       Cx - capacitance per unit length [F/m]
%       len - transmission line length [m]
%       ws - nominal system frequency [rad/s]
%       w - vector of desired frequency range
% 
% Output:
%       sys - frequency-response data model
% 
% 
% **** Comments: ****
% Inputs and outputs in are related as I(s)=Y(s)V(s), where
% I(s) = [I_{d,S}, I_{d,R}, I_{q,S}, I_{q,R}, I_{0,S}, I_{0,R}],
% V(s) = [V_{d,S}, V_{d,R}, V_{q,S}, V_{q,R}, V_{0,S}, V_{0,R}],
% and Y(s) is the nodal admittance matrix defined as 
% Y(s) = [Y'(s)/2+1/Z'(s) -1/Z'(s);
%         -1/Z'(s) Y'(s)/2+1/Z'(s)]


narginchk(6,6)

% Pre-calculation steps
wl = length(w);
CxsP = Cx*(1i.*w+1i*ws);
CxsN = Cx*(1i.*w-1i*ws);
RxLxP = Rx+Lx*(1i.*w+1i*ws);
RxLxN = Rx+Lx*(1i.*w-1i*ws);
SqrtP = sqrt(CxsP.*RxLxP);
SqrtN = sqrt(CxsN.*RxLxN);
CschP = csch(len*SqrtP);
CschN = csch(len*SqrtN);
TanhP = tanh(1/2*len*sqrt(CxsP.*RxLxP));
TanhN = tanh(1/2*len*sqrt(CxsN.*RxLxN));
denP = sqrt(RxLxP./CxsP);
denN = sqrt(RxLxN./CxsN);
arg0 = len*sqrt(Cx*1i.*w.*(Rx+Lx*1i.*w));
den0 = sqrt((Rx+Lx*1i.*w)./(Cx*1i.*w));

% Full admittance matrix
n11 = [(CschP + TanhP)./denP -CschP./denP;
    -CschP./denP (CschP + TanhP)./denP];
n12 = [(CschN + TanhN)./denN -CschN./denN;
    -CschN./denN (CschN + TanhN)./denN];
N1 = 1/2*(n11 + n12);
jN2 = 1i.*(1/2*(n11 - n12));
Ybus = [(csch(arg0) + tanh(1/2*arg0))./den0 -(csch(arg0))./den0;
    -(csch(arg0))./den0 (csch(arg0)+tanh(1/2*arg0))./den0];

% Frequency response (dq0 model)
ltlMdl = [N1 jN2 zeros(2,2*wl);
    -jN2 N1 zeros(2,2*wl);
    zeros(2,2*wl) zeros(2,2*wl) Ybus];

% Frequency responce model
Hresp = zeros(6,6,wl);
for i = 1:6
    for k = 1:6
        Hresp(i,k,:) = ltlMdl(i,((k-1)*length(w)+1):(k*wl));
    end
end
sys = frd(Hresp,w);

end
