function [cey] = energyVelTransverse(dat)
% ENERGYVELTRANSVERSE - Compute the energy velocity transverse to the waveguide.
% The transverse energy velocity is the ratio of total transverse power flux 
% to the total stored energy. The transverse energy velocity will usually be 
% nonzero for propagation in nonprincipal directions of an anisotropic material.
% 
% Usage: 
% > cey = energyVelTransverse(dat); 
% 
% Arguments: 
% - dat:    GEWdat object that stores the computed modes. 
%
% See also: ENERGYVELAXIAL, POWERFLUXTRANSVERSE.
% 
% 2022-2024 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France

if ~isscalar(dat) % compute recursively for every waveguide problem in the vector "dat"
    cey = arrayfun(@energyVelTransverse,dat,'UniformOutput',false); % apply to every object in the array "dat"
    return;
end

Py = powerFluxTransverse(dat); % orthogonal to k and plate's normal
H = energyTotal(dat);
cey = Py./H;

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
