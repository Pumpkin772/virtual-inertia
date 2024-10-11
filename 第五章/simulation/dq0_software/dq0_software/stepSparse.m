function [y, t] = stepSparse(A, B, C, D, Tf, T, display)
% STEPSPARSE Computes the step response of a state-space model.
%
% Usage:
%        stepSparse(A, B, C, D, Tf, T)
%           when invoked with no output arguments, this function plots the
%           step response on the screen
%       [y, t] = stepSparse(A, B, C, D, Tf, T)
%           returns outputs, assuming that a unit step is
%           simultaneously applied to each of the inputs
%       [y, t] = stepSparse(A, B(:,m), C, D(:,m), Tf, T)
%           returns outputs, assuming that a unit step is
%           applied to the m'th input
%       [y, t] = stepSparse(A, B(:,m), C(p,:), D(p,m), Tf, T)
%           returns the p'th output, assuming that a unit step is
%           applied to the m'th input
% 
% where
%       A, B, C, D - system matrices
%       Tf - simulation final time such that t\in[0,Tf]
%       T - numeric step size (sampling time)
%       display - computation progress is printed on screen (optional)
%
% Outputs:
%       t - time vector
%       y - the step response
% 
% See also: ssNetw


narginchk(6,7)
if nargin == 6
    display = 0;
end

NN = size(A,1);
speyeNN = speye(NN);
BB2 = sum(B,2);
DD2 = sum(D,2);

% Bilinear transformation (Tustin's method)
cond_est = condest(speyeNN-(T/2)*A);
if (cond_est>1e18)
    disp('Warning: results may be inaccurate.');
    disp('Consider reducing the step size T.');
    cond_est
    beep;
end

border_of_large_system = 500; % number of states
% For NN < border_of_large_system computation uses the explicit inverse
% For NN >= border_of_large_system computation is based on LU decomposition
if (NN < border_of_large_system )
% Solver uses the explicit inverse inv(speyeNN-(T/2)*A)
% This solver typically runs faster for smaller systems.
    
    invAT = inv(speyeNN-(T/2)*A);
% If invAT is not very sparse, and if enough memory is available,
% convert invAT to a full matrix to speed-up computations
    dummy = pi;
    dd = whos('dummy'); size_of_double=dd.bytes;
    dd = whos('invAT'); membytes = dd.bytes;
    size_of_full = prod(size(invAT))*size_of_double;
    memy = memory;
    if (size_of_full < 1.2*membytes)
        invAT = full(invAT);
    elseif ((size_of_full < 2*membytes) && ...
        (size_of_full < 0.3*memy.MaxPossibleArrayBytes))
        invAT = full(invAT);
    end
    
    Q1 = invAT*(speyeNN+(T/2)*A);
    Q2 = invAT*BB2*T;
    clear invAT;
% Test the inverse. WW should be a zero matrix:
    WW = Q1 - (speyeNN+(T/2)*A*(speyeNN+Q1)); % theoretically should be zero
% Verify that WW is indeed nearly zero:
    Wmax = full(max(max(abs(WW))));
    clear WW;
    if (Wmax > 1e-5)
        disp('Warning: step response may be inaccurate.');
        disp('Consider reducing the step size T.');
% Wmax
    end
    if (Wmax > 1e-2)
        error('stepSparse:ReduceStep', ...
        'Consider reducing the step size T.');
    end
    
% Compute step response:
    len = length(0:T:Tf);
    y = zeros(size(C,1),len);
    t = zeros(1,len);
    x = Q2/2; % initial x(t=0)
    y(:,1) = DD2;
    t(1) = 0;
    for ii=2:len
        if display
            if (mod(ii+8,10)==0)
                fprintf('Computing step response %d/%d\n',ii,len);
            end
        end
        x = Q1*x + Q2; % new value for x
        y(:,ii) = C*x+DD2; % new output
        t(ii) = (ii-1)*T;
    end
else
% When NN >= border_of_large_system computation is based on the LU decomposition.
% The inverse is never computed explicitly.
% Memory requirements are relatively low.

% LU decomposition
    [L, U, P, F, R] = lu(speyeNN-(T/2)*A);
    invR = inv(R); % diagonal matrix
    Q2 = T*(F*(U\(L\(P*invR*BB2))));
    
% Compute step response:
    len = length(0:T:Tf);
    y = zeros(size(C,1),len);
    t = zeros(1,len);
    x = Q2/2; % initial x(t=0)
    y(:,1) = DD2;
    t(1) = 0;
    
    TA2 = P*invR*(speyeNN+(T/2)*A); % still sparse ...
    
    for ii=2:len     
        x = F*(U\(L\(TA2*x))) + Q2; % new value for x
        y(:,ii) = C*x + DD2; % new output
        t(ii) = (ii-1)*T;
        if display
% Report progress
            if (mod(ii+18,20)==0)
                fprintf('Computing step response %4.2f/%4.2f\n',t(ii),Tf);
            end
        end
    end
end

if nargout == 0
    set(0,'defaulttextinterpreter','latex')
    set(0,'defaultfigurecolor',[1 1 1])
    set(0,'defaultaxesfontsize',9);
    
    plot(t,y,'LineWidth',0.5);
    xlabel('Time [s]');
    ylabel('Amplitude');
    title('\bf Step Response');
    set(gca,'TickLabelInterpreter', 'latex')
    if size(C,1)>1
        for i=1:size(C,1)
            legendInfo{i} = ['$y_{' num2str(i) '}$'];
        end
        lgd = legend(legendInfo);
        set(lgd,'Interpreter','latex')
    end
end

end
