function [ce] = energyVel(wguide, dat)
%ENERGYVEL Compute the energy velocity ce.

P = powerFlux(wguide, dat);
H = energyTotal(wguide, dat);
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
