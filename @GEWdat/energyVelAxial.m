function [cex] = energyVelAxial(dat)
% ENERGYVELAXIAL - Energy velocity component ce_x along the wave vector k.
% The energy velocity vector ce is the ratio of the total power flux vector P to the 
% total stored energy H: ce = P/H . 
% The axial component is the x-component, i.e., ce_x = P_x/H. 
% 
% Usage: 
% > cex = energyVelAxial(dat); 
% 
% Arguments: 
% - dat:    GEWdat object that stores the computed modes. 
%
% See also: ENERGYVEL, ENERGYVELTRANSVERSE, POWERFLUXAXIAL, POWERFLUX
% 
% 2022-2024 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France

if ~isscalar(dat) % compute recursively for every waveguide problem in the vector "dat"
    cex = arrayfun(@energyVelAxial,dat,'UniformOutput',false); % apply to every object in the array "dat"
    return;
end

Px = powerFluxAxial(dat); % component aligned with k-vector
H = energyTotal(dat);
cex = Px./H;

% NOTE: this function depends on many of the field calculation functions:
% - GEWintegrate
% - powerFlux
%   - poyntingVec
%       - velocity
%       - stress
% - energyTotal 
%   - energyKinetic
%       - energyDensityKinetic
%           - velocity
% 
% Additionally, without equipartion of energy, it also epends on
%   - energyElastic
%       - energyDensityElastic
%           - strain
%           - stress
% 
% Accordingly, energy velocity computations can be used to test the correct
% impementation of the field computation functions.

end
