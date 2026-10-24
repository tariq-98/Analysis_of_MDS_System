function print_results(stiffness, wn_free, undamped, damped)
%PRINT_RESULTS Summary table of the identified system parameters.

fprintf('\nStatic calibration\n');
fprintf('  Spring stiffness k (endpoints)      %10.1f N/m\n', stiffness.k);
fprintf('  Spring stiffness k (least squares)  %10.1f N/m   R^2 = %.5f\n', ...
    stiffness.k_fit, stiffness.r_squared);
fprintf('  Natural frequency, free vibration   %10.2f rad/s  (%.3f Hz)\n', ...
    wn_free, wn_free / (2 * pi));

fprintf('\n%-34s %12s %12s\n', 'Forced response', 'No damper', 'Damper');
row('Resonant peak (rad/s)', undamped.w_peak, damped.w_peak, '%12.2f');
row('Damping ratio zeta', undamped.zeta, damped.zeta, '%12.4f');
row('Effective mass (kg)', undamped.mass, damped.mass, '%12.3f');
row('Critical damping c_cr (Ns/m)', undamped.c_cr, damped.c_cr, '%12.2f');
row('Damping coefficient c (Ns/m)', undamped.c, damped.c, '%12.3f');
fprintf('\n');
end

function row(name, a, b, fmt)
fprintf(['%-34s ' fmt ' ' fmt '\n'], name, a, b);
end
