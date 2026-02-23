function [Eelec] = energyElectric(dat)
% energyElectric - Compute the total dielectric energy across the waveguide cross
% section.
%
% Usage: 
% > Eelec = energyElectric(dat);
% 
% Arguments: 
% - dat:    GEWdat object that stores the computed modes. 
%

if ~isscalar(dat) % compute recursively for every waveguide problem in the vector "dat"
    Eelec = arrayfun(@energyElectric,dat,'UniformOutput',false); % apply to every object in the array "dat"
    return;
end

eelec = energyDensityElectric(dat);
Eelec = GEWintegrate(dat.gew, eelec); % integrate over all layers

end
