function [Eelast] = energyElastic(wguide,dat)
%ENERGYELASTIC Compute the total elastic energy across the waveguide cross
%section.

eelast = energyDensityElastic(wguide, dat);
Eelast = GEWintegrate(wguide, eelast); % integrate over all layers

end
