function [eElastic] = energyDensityElastic(dat)
%ENERGYDENSITYELASTIC Compute the elastic energy density.
% 
% eElastic = 1/4*real(T:S*)
%
% Usage: 
% > eElastic = energyDensityElastic(dat);
% 
% Arguments: 
% - dat:    GEWdat object that stores the computed modes. 
%

if ~isscalar(dat) % compute recursively for every waveguide problem in the vector "dat"
    eElastic = arrayfun(@energyDensityElastic,dat,'UniformOutput',false); % apply to every object in the array "dat"
    return;
end

if isDissipative(dat.gew)
    warning('energyDensityElastic(): Is it valid for viscous media?')
end

S = strain(dat);
T = stress(dat); 

eElastic = cell(dat.gew.geom.nLay,1);
for l = 1:dat.gew.geom.nLay
    eElastic{l} = 1/4*real(sum(sum(conj(S{l}).*T{l}, 5), 4)); % double contraction
end

end
