function data = load_lab_data(data_dir)
%LOAD_LAB_DATA Read the four experiment CSV files into a struct.
%
%   data = LOAD_LAB_DATA(data_dir) returns
%     data.static.load_N, data.static.sensor_mV     static calibration
%     data.periods_ms                               free-vibration periods
%     data.undamped.f_Hz, data.undamped.amp_mV      forced sweep, no damper
%     data.damped.f_Hz,   data.damped.amp_mV        forced sweep, damper fitted

static = read_csv(fullfile(data_dir, 'static_calibration.csv'));
data.static.load_N = static(:, 1).';
data.static.sensor_mV = static(:, 2).';

periods = read_csv(fullfile(data_dir, 'free_vibration_periods.csv'));
data.periods_ms = periods(:, 1).';

undamped = read_csv(fullfile(data_dir, 'sweep_undamped.csv'));
data.undamped.f_Hz = undamped(:, 1).';
data.undamped.amp_mV = undamped(:, 2).';

damped = read_csv(fullfile(data_dir, 'sweep_damped.csv'));
data.damped.f_Hz = damped(:, 1).';
data.damped.amp_mV = damped(:, 2).';
end

function m = read_csv(path)
% dlmread skips the header row and works in both MATLAB and GNU Octave.
if ~exist(path, 'file')
    error('load_lab_data:missingFile', 'Data file not found: %s', path);
end
m = dlmread(path, ',', 1, 0);
end
