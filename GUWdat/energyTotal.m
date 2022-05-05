function [Etot] = energyTotal(wguide,dat)
%ENERGYTOTAL Compute total average stored energy across waveguide cross section.
% Note: Equipartition of energy is not valid for waves with complex wavenumbers.
% For this reason the total energy is explicitly computed as the sum of elastic
% and kinetik energies. 
% see: 

Etot = energyKinetic(wguide,dat) + energyElastic(wguide,dat);

end
