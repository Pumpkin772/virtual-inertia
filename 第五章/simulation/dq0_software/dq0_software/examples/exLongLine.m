% Demonstrates how to construct analytic and approximate models of
% a long transmission line in the dq0 reference frame

% Default settings for figures
set(0,'defaulttextinterpreter','latex')
set(0,'defaultfigurecolor',[1 1 1])
set(0,'defaultaxesfontsize',9);

% Parameters from Stevenson's book p.99
mi = 1609.344; % 1mi in m
Rx = 0.172/mi; % [Ohm/m]
Lx = 2.18/mi*1e-3; % [H/m]
Cx = 0.0136/mi*1e-6; % [F/m]
len = 225*mi; % [m]
ws = 2*pi*60;
w = logspace(0,6,1000);


% Construct the analytic dq0 model
sysAnalyt = longLineAnalyt(Rx, Lx, Cx, len, ws, w);
% Construct the approximate (N=10) dq0 model
[At, Bt, Ct, Dt] = longLine(Rx, Lx, Cx, len, 10, ws);
sysApprox = ss(full(At), full(Bt), full(Ct), full(Dt));

% Compute frequency response
[mag_analyt,ph_analyt,w_analyt] = bode(sysAnalyt(1,1),w);
[mag_approx,ph_approx,w_approx] = bode(sysApprox(1,1),w);

% Draw bode plots for V_{d,S} to I_{d,S}
figure(1);
subplot(2,1,1);
semilogx(w_analyt,20*log10(squeeze(mag_analyt(1,1,:))),'k','LineWidth',0.6,...
    'Color',[0.7 0.7 0.7]);
hold on;
semilogx(w_approx,20*log10(squeeze(mag_approx(1,1,:))),'--k','LineWidth',0.5,...
    'Color',[0 0 0]);
title('$V_{d,S}\to I_{d,S}$','FontSize',9)
ylabel('Magnitude (dB)','FontSize',9)
xlim([1 1e5])
lgd = legend('analytic', 'approximate', 'Location', 'NorthWest');
set(lgd,'Interpreter','latex', 'FontSize',9)

subplot(2,1,2);
semilogx(w_analyt,squeeze(ph_analyt(1,1,:)),'k','LineWidth',0.6,...
    'Color',[0.7 0.7 0.7]);
hold on;
semilogx(w_approx,squeeze(ph_approx(1,1,:)),'--k','LineWidth',0.5,...
    'Color',[0 0 0]);
xlabel('Frequency [rad/s]','FontSize',9)
ylabel('Phase (deg)','FontSize',9)
xlim([1 1e5])

axesHandles = findall(0,'type','axes');
set(axesHandles,'TickLabelInterpreter', 'latex',...
    'XTick', [10.^0 10.^3 10.^4 10.^5],'xgrid','on','ygrid','on')
