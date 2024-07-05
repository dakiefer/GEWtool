function [Eelast] = energyElastic(gew,dat)
%ENERGYELASTIC Compute the total elastic energy across the waveguide cross
%section.

if ~isscalar(gew) % compute recursively for every waveguide problem in the vector "gew"
    compute = @(gewObj,datObj) energyElastic(gewObj, datObj); % function to apply
    Eelast = arrayfun(compute,gew,dat,'UniformOutput',false); % apply to every object in the arrays "gew" and "dat"
    return;
end

eelast = energyDensityElastic(gew, dat);
Eelast = GEWintegrate(gew, eelast); % integrate over all layers

end
