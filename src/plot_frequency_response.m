function fig = plot_frequency_response(omega, amplitude, hp, label)
%PLOT_FREQUENCY_RESPONSE Swept amplitude with the half-power bandwidth marked.

fig = figure('Name', ['Frequency response, ' label]);
plot(omega, amplitude, '-ok', 'LineWidth', 1.5, 'MarkerSize', 4);
hold on
plot([hp.w1 hp.w2], [hp.half_power hp.half_power], '--r', 'LineWidth', 1.2);
plot([hp.w1 hp.w1], [0 hp.half_power], ':r');
plot([hp.w2 hp.w2], [0 hp.half_power], ':r');
plot(hp.w_peak, hp.amp_peak, 'xr', 'MarkerSize', 10, 'LineWidth', 2);
hold off
grid on
xlabel('Forcing frequency \omega (rad/s)');
ylabel('Amplitude (mV)');
title(sprintf('Frequency response, %s: \\zeta = %.4f', label, hp.zeta));
legend('Measured', 'Half-power bandwidth', 'Location', 'northwest');
end
