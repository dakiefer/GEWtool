function [P] = powerFlux(dat)
% POWERFLUX - Power flux vectors.
%
% Computes the total power flux vector P. This is done by integrating the power
% flux density p (poynting vector) over the waveguide's cross-section (e.g., plate
% thickness).
% 
% Usage: 
% > P = powerFlux(dat);
% 
% Arguments: 
% - dat:    GEWdat object that stores the computed modes. 
%
% See also: powerFluxAxial, powerFluxTransverse
% 
% 2022-2024 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France

if ~isscalar(dat) % compute recursively for every waveguide problem in the vector "dat"
    P = arrayfun(@powerFlux,dat,'UniformOutput',false); % apply to every object in the array "dat"
    return;
end

if isa(dat.gew,"CylinderCircumferential")
    warning('Circumferential waves do not support this function yet. The results will be wrong.');
end
p = poyntingVec(dat);
P = GEWintegrate(dat.gew, p);

end
