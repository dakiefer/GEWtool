function [Ekin] = energyKinetic(wguide, dat)
%ENERGYKIN Compute the total kinetik energy over the waveguide cross section.

ekin = energyDensityKinetic(wguide, dat);
Ekin = GEWintegrate(wguide, ekin);

end
