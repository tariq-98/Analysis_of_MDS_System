function fig = plot_phase_lag(omega, phi, wn)
%PLOT_PHASE_LAG Phase lag of the damped system against forcing frequency.
%   The phase is not measured on the rig: it is evaluated from the SDOF
%   model at the sweep frequencies using the identified zeta and wn.

fig = figure('Name', 'Phase lag');
plot(omega, phi, '-ok', 'LineWidth', 1.5, 'MarkerSize', 4);
hold on
plot([wn wn], [-pi 0], ':r');
plot([min(omega) max(omega)], [-pi/2 -pi/2], ':r');
hold off
grid on
ylim([-pi 0]);
set(gca, 'YTick', [-pi -pi/2 0], 'YTickLabel', {'-\pi', '-\pi/2', '0'});
xlabel('Forcing frequency \omega (rad/s)');
ylabel('Phase lag \phi (rad)');
title('Predicted phase lag with damper fitted (from identified \zeta, \omega_n)');
end
