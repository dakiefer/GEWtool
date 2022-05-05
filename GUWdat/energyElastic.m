function [Eelast] = energyElastic(wguide,dat)
%ENERGYELASTIC Compute the total elastic energy across the waveguide cross
%section.

eelast = energyDensityElastic(wguide, dat);

Eelast = zeros(size(dat.k));
for i = 1:wguide.geom.nLay % for every layer
    lims = wguide.geom.yItf(i,:);   % lower and upper integration limits
    Eelast = Eelast + chebintegrate(eelast{i}, lims, 3);
end

end
