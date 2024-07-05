function [eElastic] = energyDensityElastic(gew, dat)
%ENERGYDENSITYELASTIC Compute the elastic energy density.

if ~isscalar(gew) % compute recursively for every waveguide problem in the vector "gew"
    compute = @(gewObj,datObj) energyDensityElastic(gewObj, datObj); % function to apply
    eElastic = arrayfun(compute,gew,dat,'UniformOutput',false); % apply to every object in the arrays "gew" and "dat"
    return;
end

if isDissipative(gew)
    warning('energyDensityElastic(): Is it valid for viscous media?')
end

S = strain(gew, dat);
T = stress(gew, dat);

eElastic = cell(gew.geom.nLay,1);
for i = 1:gew.geom.nLay
    eElastic{i} = 1/4*real(sum(sum(conj(S{i}).*T{i}, 5), 4)); % double contraction
end

end
