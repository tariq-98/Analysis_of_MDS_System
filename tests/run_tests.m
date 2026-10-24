function run_tests()
%RUN_TESTS Test suite for the MDS analysis. Runs in MATLAB and GNU Octave.
%
%   From the repository root:   addpath('tests'); run_tests
%   Exits with an error (non-zero status under `octave --eval`) if any fail.

root = fileparts(fileparts(mfilename('fullpath')));
addpath(root, fullfile(root, 'src'));

tests = {
    @natural_frequency_from_known_period
    @stiffness_recovers_known_spring
    @half_power_recovers_damping_of_analytic_response
    @half_power_is_independent_of_noise_away_from_resonance
    @half_power_rejects_unbracketed_sweep
    @half_power_rejects_mismatched_inputs
    @phase_lag_passes_through_minus_half_pi_at_resonance
    @system_properties_follow_definitions
    @pipeline_reproduces_original_lab_results
};

failures = 0;
for i = 1:numel(tests)
    name = func2str(tests{i});
    try
        tests{i}();
        fprintf('  pass  %s\n', name);
    catch err
        failures = failures + 1;
        fprintf('  FAIL  %s\n        %s\n', name, err.message);
    end
end
fprintf('\n%d of %d tests passed\n', numel(tests) - failures, numel(tests));
if failures > 0
    error('run_tests:failed', '%d test(s) failed', failures);
end
end

% ---------------------------------------------------------------- tests

function natural_frequency_from_known_period()
assert_close(natural_frequency([100 100 100]), 2 * pi / 0.1, 1e-12);
end

function stiffness_recovers_known_spring()
rig = rig_constants();
k_true = 2500;
load_N = 0:5:25;
extension_mm = load_N / k_true * 1e3;
sensor_mV = extension_mm * rig.AB_mm / rig.AC_mm * rig.sensor_mV_per_mm;
s = spring_stiffness(load_N, sensor_mV, rig);
assert_close(s.k, k_true, 1e-9);
assert_close(s.k_fit, k_true, 1e-9);
assert_close(s.r_squared, 1, 1e-12);
end

function half_power_recovers_damping_of_analytic_response()
% Densely sampled SDOF response with a known damping ratio.
wn = 45;
for zeta = [0.01 0.02 0.05]
    omega = linspace(0.5 * wn, 1.5 * wn, 4001);
    hp = half_power_bandwidth(omega, sdof_amplitude(omega, wn, zeta));
    assert_close(hp.zeta, zeta, 0.02 * zeta);
    assert_close(hp.w_peak, wn * sqrt(1 - 2 * zeta ^ 2), 0.01);
end
end

function half_power_is_independent_of_noise_away_from_resonance()
% Bumps far from the peak must not move the result: only the points that
% bracket the half-power crossings matter.
wn = 45; zeta = 0.02;
omega = linspace(0.5 * wn, 1.5 * wn, 2001);
clean = sdof_amplitude(omega, wn, zeta);
noisy = clean;
far = abs(omega - wn) > 0.2 * wn;
noisy(far) = noisy(far) .* (1 + 0.2 * sin(37 * omega(far)));
assert_close(half_power_bandwidth(omega, noisy).zeta, ...
    half_power_bandwidth(omega, clean).zeta, 1e-12);
end

function half_power_rejects_unbracketed_sweep()
omega = linspace(40, 46, 50);           % stops at the peak
assert_error(@() half_power_bandwidth(omega, sdof_amplitude(omega, 46, 0.02)), ...
    'half_power_bandwidth:notBracketed');
end

function half_power_rejects_mismatched_inputs()
assert_error(@() half_power_bandwidth(1:5, 1:4), 'half_power_bandwidth:sizeMismatch');
end

function phase_lag_passes_through_minus_half_pi_at_resonance()
wn = 43; zeta = 0.05;
assert_close(phase_lag(wn, wn, zeta), -pi / 2, 1e-12);
assert_close(phase_lag(1e-6, wn, zeta), 0, 1e-6);
assert_close(phase_lag(1e9, wn, zeta), -pi, 1e-6);
phi = phase_lag(linspace(1, 100, 500), wn, zeta);
assert(all(diff(phi) < 0), 'phase lag should fall monotonically');
end

function system_properties_follow_definitions()
s = system_properties(3600, 60, 0.1);
assert_close(s.mass, 1, 1e-12);
assert_close(s.c_cr, 120, 1e-12);
assert_close(s.c, 12, 1e-12);
end

function pipeline_reproduces_original_lab_results()
% Reference values produced by the original single-script analysis
% (Experiment1.m, 2020), run under GNU Octave 8.4 before the refactor.
r = run_analysis('Plot', false, 'Verbose', false);
tol = 1e-9;
assert_rel(r.stiffness.k, 3411.114542, tol);
assert_rel(r.wn_free, 45.56334523, tol);
assert_rel(r.undamped.w_peak, 47.1238898, tol);
assert_rel(r.undamped.zeta, 0.02930940475, tol);
assert_rel(r.undamped.mass, 1.643104172, tol);
assert_rel(r.undamped.c_cr, 149.7306453, tol);
assert_rel(r.undamped.c, 4.388516087, tol);
assert_rel(r.damped.w_peak, 42.72566009, tol);
assert_rel(r.damped.zeta, 0.01363505458, tol);
assert_rel(r.damped.mass, 1.86861031, tol);
assert_rel(r.damped.c_cr, 159.6752179, tol);
assert_rel(r.damped.c, 2.177180311, tol);
expected_phase = [-0.0115922 -0.0149374 -0.0190924 -0.0245238 -0.0320961 ...
    -0.043625 -0.0636975 -0.10823 -0.293373 -0.743927 -1.5708 -2.79312 ...
    -2.90754 -3.00354 -3.03261 -3.05826 -3.08107 -3.09877];
assert_close(r.phase_lag, expected_phase, 1e-5);
end

% -------------------------------------------------------------- helpers

function a = sdof_amplitude(omega, wn, zeta)
r = omega / wn;
a = 1 ./ sqrt((1 - r .^ 2) .^ 2 + (2 * zeta * r) .^ 2);
end

function assert_close(actual, expected, tol)
err = max(abs(actual(:) - expected(:)));
if ~(err <= tol)
    error('expected %s, got %s (error %g > %g)', ...
        mat2str(expected, 8), mat2str(actual, 8), err, tol);
end
end

function assert_rel(actual, expected, tol)
assert_close(actual, expected, tol * abs(expected));
end

function assert_error(fn, id)
try
    fn();
catch err
    if ~strcmp(err.identifier, id)
        error('expected error %s, got %s', id, err.identifier);
    end
    return
end
error('expected error %s, but none was raised', id);
end
