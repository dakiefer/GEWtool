function [ce] = energyVel(gew, dat)
%ENERGYVEL Compute the energy velocity ce.

P = powerFlux(gew, dat);
H = energyTotal(gew, dat);
ce = P./H;

% NOTE: this function depends on basically all field calculation functions:
% - GEWintegrate
% - powerFlux
%   - poyntingVec
%       - velocity
%       - stress
% - energyTotal
%   - energyKinetic
%       - energyDensityKinetic
%           - velocity
%   - energyElastic
%       - energyDensityElastic
%           - strain
%           - stress

end
