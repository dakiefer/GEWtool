function [Ekin] = energyKinetic(wguide, dat)
%ENERGYKIN Compute the total kinetik energy over the waveguide cross section.

ekin = energyDensityKinetic(wguide, dat);

Ekin = zeros(size(dat.k));
for i = 1:wguide.geom.nLay
    lims = wguide.geom.yItf(i,:);   % lower and upper integration limits
    Ekin = Ekin + chebintegrate(ekin{i}, lims, 3);
end

end
