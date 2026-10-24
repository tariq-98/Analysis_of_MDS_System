function result = spring_stiffness(load_N, sensor_mV, rig)
%SPRING_STIFFNESS Spring constant from the static load test.
%
%   result = SPRING_STIFFNESS(load_N, sensor_mV, rig) converts transducer
%   readings to spring extension through the lever geometry, then estimates
%   the stiffness two ways:
%
%     result.k           endpoint slope, first to last reading (N/m). This is
%                        the value used through the rest of the analysis.
%     result.k_fit       least-squares slope through the origin (N/m)
%     result.r_squared   coefficient of determination of that fit
%     result.extension_m spring extension at each load (m)

travel_B_mm = sensor_mV / rig.sensor_mV_per_mm;
extension_mm = travel_B_mm * rig.AC_mm / rig.AB_mm;
extension_m = extension_mm * 1e-3;

result.extension_m = extension_m;
result.k = (load_N(end) - load_N(1)) / (extension_m(end) - extension_m(1));

% Hooke's law passes through the origin, so fit F = k*x with no intercept.
x = extension_m(:);
F = load_N(:);
result.k_fit = (x' * F) / (x' * x);
residual = F - result.k_fit * x;
result.r_squared = 1 - sum(residual .^ 2) / sum((F - mean(F)) .^ 2);
end
