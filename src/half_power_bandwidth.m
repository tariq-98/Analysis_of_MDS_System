function hp = half_power_bandwidth(omega, amplitude)
%HALF_POWER_BANDWIDTH Damping ratio from a frequency response by the
%half-power (-3 dB) method.
%
%   hp = HALF_POWER_BANDWIDTH(omega, amplitude) takes a swept response,
%   finds the resonant peak, and walks outwards on each side until the
%   amplitude falls below peak/sqrt(2). The crossing frequencies w1 and w2
%   are found by linear interpolation between the bracketing points, and
%
%       zeta = (w2 - w1) / (2 * w_peak)
%
%   Returns a struct with fields
%     w_peak, amp_peak   resonant frequency (rad/s) and amplitude
%     half_power         amplitude threshold, amp_peak / sqrt(2)
%     w1, w2             lower and upper half-power frequencies (rad/s)
%     zeta               damping ratio
%
%   Walking out from the peak, rather than interpolating over a fixed index
%   window, means the method copes with noisy, non-monotonic data away
%   from resonance and needs no hand-picked indices per data set.

omega = omega(:).';
amplitude = amplitude(:).';
if numel(omega) ~= numel(amplitude)
    error('half_power_bandwidth:sizeMismatch', ...
        'omega has %d points but amplitude has %d', numel(omega), numel(amplitude));
end

[amp_peak, ip] = max(amplitude);
threshold = amp_peak / sqrt(2);

% Lower side: last point below the threshold before the peak.
il = find(amplitude(1:ip) < threshold, 1, 'last');
% Upper side: first point below the threshold after the peak.
iu = ip - 1 + find(amplitude(ip:end) < threshold, 1, 'first');
if isempty(il) || isempty(iu)
    error('half_power_bandwidth:notBracketed', ...
        'The sweep does not fall below the half-power level on both sides of the peak.');
end

hp.w_peak = omega(ip);
hp.amp_peak = amp_peak;
hp.half_power = threshold;
hp.w1 = crossing(omega(il), amplitude(il), omega(il + 1), amplitude(il + 1), threshold);
hp.w2 = crossing(omega(iu - 1), amplitude(iu - 1), omega(iu), amplitude(iu), threshold);
hp.zeta = (hp.w2 - hp.w1) / (2 * hp.w_peak);
end

function w = crossing(wa, ya, wb, yb, level)
w = wa + (level - ya) * (wb - wa) / (yb - ya);
end
