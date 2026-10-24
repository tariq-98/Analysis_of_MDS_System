function rig = rig_constants()
%RIG_CONSTANTS Geometry and sensor calibration of the vibration test rig.
%
%   The displacement transducer sits at point B on a pivoted beam, while the
%   spring acts at point C. Spring extension is recovered from the sensor
%   reading by the lever ratio AC/AB.
%
%   rig.sensor_mV_per_mm  transducer sensitivity (mV per mm of travel at B)
%   rig.AB_mm             pivot to transducer distance
%   rig.AC_mm             pivot to spring distance

rig.sensor_mV_per_mm = 350;
rig.AB_mm = 125;
rig.AC_mm = 690;
end
