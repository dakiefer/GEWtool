function [Eelast] = energyElastic(gew,dat)
%ENERGYELASTIC Compute the total elastic energy across the waveguide cross
%section.

eelast = energyDensityElastic(gew, dat);
Eelast = GEWintegrate(gew, eelast); % integrate over all layers

end
