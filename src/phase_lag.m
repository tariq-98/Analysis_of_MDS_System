function phi = phase_lag(omega, wn, zeta)
%PHASE_LAG Phase of a damped SDOF response relative to the forcing (rad).
%
%   phi = PHASE_LAG(omega, wn, zeta) evaluates
%
%       phi = -atan2(2*zeta*r, 1 - r^2),   r = omega / wn
%
%   which runs from 0 well below resonance, through -pi/2 at resonance,
%   towards -pi well above it. atan2 keeps the result on the correct
%   branch without a manual quadrant correction.

r = omega / wn;
phi = -atan2(2 * zeta * r, 1 - r .^ 2);
end
