function [Eelast] = energyElastic(wguide,dat)
%ENERGYELASTIC Compute the total elastic energy across the waveguide cross
%section.

eelast = energyDensityElastic(wguide, dat);
Eelast = chebintegrate(eelast, wguide.geom.yItf, 3);

end
