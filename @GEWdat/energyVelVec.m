function [ce] = energyVelVec(dat)
% ENERGYVEL - Energy velocity vector ce.
% The energy velocity vector ce is the ratio of the total power flux vector P to the 
% total stored energy H: ce = P/H . 
% 
% Usage: 
% > ce = energyVelVec(dat); 
% 
% Arguments: 
% - dat:    GEWdat object that stores the computed modes. 
%
% See also: energyVelAxial, energyVelTransverse, powerFlux
% 
% 2024 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France

if ~isscalar(dat) % compute recursively for every waveguide problem in the vector "dat"
    ce = arrayfun(@energyVelVec,dat,'UniformOutput',false); % apply to every object in the array "dat"
    return;
end

P = powerFlux(dat);   % for a plate: vector in the plane of the plate
H = energyTotal(dat);
ce = P./H*dat.gew.np.fh0;

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
