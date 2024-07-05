function [Ekin] = energyKinetic(gew, dat)
%ENERGYKIN Compute the total kinetik energy over the waveguide cross section.

if ~isscalar(gew) % compute recursively for every waveguide problem in the vector "gew"
    compute = @(gewObj,datObj) energyKinetic(gewObj, datObj); % function to apply
    Ekin = arrayfun(compute,gew,dat,'UniformOutput',false); % apply to every object in the arrays "gew" and "dat"
    return;
end

ekin = energyDensityKinetic(gew, dat);
Ekin = GEWintegrate(gew, ekin);

end
