function [Etot] = energyTotal(wguide,dat)
%ENERGYTOTAL Compute total average stored energy across waveguide cross section.

Etot = energyKinetic(wguide,dat) + energyElastic(wguide,dat);

end
