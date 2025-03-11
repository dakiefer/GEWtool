function [Ekin] = energyKinetic(dat)
%ENERGYKIN Compute the total kinetik energy over the waveguide cross section.
%
% Usage: 
% > Ekin = energyKinetic(dat);
% 
% Arguments: 
% - dat:    GEWdat object that stores the computed modes. 
%

if ~isscalar(dat) % compute recursively for every waveguide problem in the vector "dat"
    Ekin = arrayfun(@energyKinetic,dat,'UniformOutput',false); % apply to every object in the array "dat"
    return;
end

ekin = energyDensityKinetic(dat);
Ekin = GEWintegrate(dat.gew, ekin);

end
