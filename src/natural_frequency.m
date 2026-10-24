function wn = natural_frequency(periods_ms)
%NATURAL_FREQUENCY Natural frequency (rad/s) from free-vibration periods.
%
%   wn = NATURAL_FREQUENCY(periods_ms) averages repeated period measurements
%   in milliseconds and returns wn = 2*pi / T.

T = mean(periods_ms) * 1e-3;
wn = 2 * pi / T;
end
