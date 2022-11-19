function [Ekin] = energyKinetic(gew, dat)
%ENERGYKIN Compute the total kinetik energy over the waveguide cross section.

ekin = energyDensityKinetic(gew, dat);
Ekin = GEWintegrate(gew, ekin);

end
