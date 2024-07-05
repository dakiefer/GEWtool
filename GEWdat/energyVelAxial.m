function [cex] = energyVelAxial(gew, dat)
% ENERGYVELAXIAL - Energy velocity component ce_x along the wave vector k.
% The energy velocity vector ce is the radio of the total power flux vector P to the 
% total stored energy H: ce = P/H . 
% The axial component is the x-component, i.e., ce_x = P_x/H. 
% 
% See also: ENERGYVEL, ENERGYVELTRANSVERSE, POWERFLUXAXIAL, POWERFLUX
% 
% 2022-2024 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France

if ~isscalar(gew) % compute recursively for every waveguide problem in the vector "gew"
    compute = @(gewObj,datObj) energyVelAxial(gewObj, datObj); % function to apply
    cex = arrayfun(compute,gew,dat,'UniformOutput',false); % apply to every object in the arrays "gew" and "dat"
    return;
end

Px = powerFluxAxial(gew, dat); % component aligned with k-vector
H = energyTotal(gew, dat);
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
