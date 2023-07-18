function [cez] = energyVelTransverse(gew, dat)
% ENERGYVELTRANSVERSE - Compute the energy velocity transverse to the waveguide.
% The transverse energy velocity is the radio of total transverse power flux 
% to the total stored energy. The transverse energy velocity will usually be 
% nonzero for propagation in nonprincipal directions of an anisotropic material.
% 
% See also: ENERGYVELAXIAL, POWERFLUXTRANSVERSE.
% 
% 2022-2023 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France

Pz = powerFluxTransverse(gew, dat);
H = energyTotal(gew, dat);
cez = Pz./H;

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
