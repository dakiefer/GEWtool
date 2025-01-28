function [Eelast] = energyElastic(dat)
%ENERGYELASTIC Compute the total elastic energy across the waveguide cross
%section.
%
% Usage: 
% > Eelast = energyElastic(dat);
% 
% Arguments: 
% - dat:    GEWdat object that stores the computed modes. 
%

if ~isscalar(dat) % compute recursively for every waveguide problem in the vector "dat"
    Eelast = arrayfun(@energyElastic,dat,'UniformOutput',false); % apply to every object in the arrays "dat"
    return;
end

eelast = energyDensityElastic(dat);
Eelast = GEWintegrate(dat.gew, eelast); % integrate over all layers

end
