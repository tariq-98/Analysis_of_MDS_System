function results = run_analysis(varargin)
%RUN_ANALYSIS Identify the parameters of the mass-damper-spring rig.
%
%   results = RUN_ANALYSIS() loads the measurements in data/, runs the full
%   identification, prints a summary and draws the four figures.
%
%   RUN_ANALYSIS('SaveFigures', true) also writes the figures as PNGs to
%   docs/figures/.  RUN_ANALYSIS('Plot', false) skips plotting, which is
%   what the tests use.
%
%   Pipeline
%     1. Static load test        -> spring stiffness k
%     2. Free vibration periods  -> natural frequency wn
%     3. Forced sweep, no damper -> damping ratio by the half-power method
%     4. Forced sweep, damper    -> damping ratio, phase lag
%     5. k, wn, zeta             -> mass, critical and actual damping

opts = parse_options(varargin{:});
root = fileparts(mfilename('fullpath'));
addpath(fullfile(root, 'src'));

data = load_lab_data(fullfile(root, 'data'));
rig = rig_constants();

% 1-2. Static and free-vibration tests.
stiffness = spring_stiffness(data.static.load_N, data.static.sensor_mV, rig);
k = stiffness.k;
wn_free = natural_frequency(data.periods_ms);

% 3. Forced response without the damper. Mass comes from the free-vibration
%    natural frequency, which is not distorted by forcing.
omega_u = 2 * pi * data.undamped.f_Hz;
hp_u = half_power_bandwidth(omega_u, data.undamped.amp_mV);
undamped = system_properties(k, wn_free, hp_u.zeta);
undamped.w_peak = hp_u.w_peak;

% 4. Forced response with the damper. Fitting the damper adds moving parts,
%    so the natural frequency is taken from this sweep's own resonant peak
%    rather than the free-vibration value.
omega_d = 2 * pi * data.damped.f_Hz;
hp_d = half_power_bandwidth(omega_d, data.damped.amp_mV);
damped = system_properties(k, hp_d.w_peak, hp_d.zeta);
damped.w_peak = hp_d.w_peak;
phi = phase_lag(omega_d, hp_d.w_peak, hp_d.zeta);

results = struct('stiffness', stiffness, 'wn_free', wn_free, ...
    'undamped', undamped, 'damped', damped, ...
    'half_power_undamped', hp_u, 'half_power_damped', hp_d, ...
    'phase_lag', phi, 'omega_damped', omega_d);

if opts.Verbose
    print_results(stiffness, wn_free, undamped, damped);
end

if opts.Plot
    figs = [plot_load_extension(data.static.load_N, stiffness), ...
            plot_frequency_response(omega_u, data.undamped.amp_mV, hp_u, 'no damper'), ...
            plot_frequency_response(omega_d, data.damped.amp_mV, hp_d, 'with damper'), ...
            plot_phase_lag(omega_d, phi, hp_d.w_peak)];
    if opts.SaveFigures
        names = {'load_extension', 'response_undamped', 'response_damped', 'phase_lag'};
        out_dir = fullfile(root, 'docs', 'figures');
        if ~exist(out_dir, 'dir'), mkdir(out_dir); end
        for i = 1:numel(figs)
            print(figs(i), fullfile(out_dir, [names{i} '.png']), '-dpng', '-r150');
        end
    end
end
end

function opts = parse_options(varargin)
opts = struct('Plot', true, 'SaveFigures', false, 'Verbose', true);
if mod(numel(varargin), 2) ~= 0
    error('run_analysis:badOptions', 'Options must be name-value pairs.');
end
for i = 1:2:numel(varargin)
    name = varargin{i};
    if ~isfield(opts, name)
        error('run_analysis:badOptions', 'Unknown option ''%s''.', name);
    end
    opts.(name) = logical(varargin{i + 1});
end
end
