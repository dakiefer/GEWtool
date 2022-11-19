function [eElastic] = energyDensityElastic(gew, dat)
%ENERGYDENSITYELASTIC Compute the elastic energy density.

warning('energyDensityElastic(): Is it valid for viscous media?')

S = strain(gew, dat);
T = stress(gew, dat);

eElastic = cell(gew.geom.nLay,1);
for i = 1:gew.geom.nLay
    eElastic{i} = 1/4*real(sum(sum(conj(S{i}).*T{i}, 5), 4)); % double contraction
end

end
