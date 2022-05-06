function [eElastic] = energyDensityElastic(wguide, dat)
%ENERGYDENSITYELASTIC Compute the elastic energy density.

warning('energyDensityElastic(): not valid for viscous media?')

S = strain(wguide, dat);
T = stress(wguide, dat);

eElastic = cell(wguide.geom.nLay,1);
for i = 1:wguide.geom.nLay
    eElastic{i} = 1/4*real(sum(sum(conj(S{i}).*T{i}, 5), 4)); % double contraction
end

end
