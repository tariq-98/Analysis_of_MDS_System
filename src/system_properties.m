function sys = system_properties(k, wn, zeta)
%SYSTEM_PROPERTIES Mass and damping of a single-degree-of-freedom system.
%
%   sys = SYSTEM_PROPERTIES(k, wn, zeta) uses wn = sqrt(k/m) to recover
%   the effective mass, then
%
%       c_cr = 2*sqrt(m*k)      critical damping coefficient (Ns/m)
%       c    = zeta * c_cr      damping coefficient (Ns/m)

sys.k = k;
sys.wn = wn;
sys.zeta = zeta;
sys.mass = k / wn ^ 2;
sys.c_cr = 2 * sqrt(sys.mass * k);
sys.c = zeta * sys.c_cr;
end
