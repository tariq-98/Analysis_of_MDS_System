function fig = plot_load_extension(load_N, stiffness)
%PLOT_LOAD_EXTENSION Static load against spring extension, with the fit.

fig = figure('Name', 'Load vs extension');
x = stiffness.extension_m;
plot(x * 1e3, load_N, 'ok', 'MarkerFaceColor', 'k');
hold on
xs = linspace(0, max(x), 50);
plot(xs * 1e3, stiffness.k_fit * xs, '-b', 'LineWidth', 1.5);
hold off
grid on
xlabel('Spring extension (mm)');
ylabel('Load (N)');
title(sprintf('Static calibration: k = %.0f N/m (R^2 = %.4f)', ...
    stiffness.k_fit, stiffness.r_squared));
legend('Measured', 'Least-squares fit', 'Location', 'northwest');
end
