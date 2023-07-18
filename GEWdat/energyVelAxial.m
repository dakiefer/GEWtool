function [cex] = energyVelAxial(gew, dat)
% ENERGYVELAXIAL - Compute the energy velocity ce along the waveguide.
% The energy velocity is the radio of total axial power flux to the total 
% stored energy. 
% 
% See also: ENERGYVELTRANSVERSE, POWERFLUXAXIAL.
% 
% 2022-2023 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France

P = powerFlux(gew, dat);
H = energyTotal(gew, dat);
cex = P./H;

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
