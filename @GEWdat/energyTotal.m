function [Etot] = energyTotal(dat)
% ENERGYTOTAL - Total time-averaged stored energy across waveguide cross section.
% The total energy H is the sum of the kinetic energy K and the elastic energy E: 
% H = K + E .
% When the waveguide is nondissipative and only real wavenumbers are considered,
% there is equipartition of energy [1]. Then, H = 2*K holds and yields faster
% computations as we avoid computing the strain/stress tensors. 
%
% Note: Equipartition of energy is not valid for waves with complex wavenumbers.
%
% Literature: 
% [1] K.-J. Langenberg, R. Marklein, and K. Mayer, Ultrasonic Nondestructive
% Testing of Materials: Theoretical Foundations (translated from German), 1st
% ed. Boca Raton: CRC Press, 2012. doi: 10.1201/b11724.
% 
% Usage: 
% > Etot = energyTotal(dat);
% 
% Arguments: 
% - dat:    GEWdat object that stores the computed modes. 
%
% 2024 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France

if ~isscalar(dat) % compute recursively for every waveguide problem in the vector "dat"
    compute = @(gewObj,datObj) energyTotal(gewObj, datObj); % function to apply
    Etot = arrayfun(@energyTotal,dat,'UniformOutput',false); % apply to every object in the arrays "dat"
    return;
end

if ~dat.gew.isDissipative && isreal(dat.k)
    Etot = 2*energyKinetic(dat); % exploit equipartition of energy
else
    Etot = energyKinetic(dat) + energyElastic(dat);
end

end
