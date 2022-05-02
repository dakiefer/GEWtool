function [eElastic] = energyDensityElastic(wguide, dat)
%ENERGYDENSITYELASTIC Compute the elastic energy density.

warning('energyDensityElastic(): not valid for viscuous media.')

S = strain(wguide, dat);
T = stress(wguide, dat);
eElastic = 1/4*real(sum(sum(conj(S).*T, 5), 4));

end
