function [P] = powerFlux(gew, dat)
% POWERFLUX - Power flux vectors.
%
% Computes the total power flux vector P. This is done by integrating the power
% flux density p (poynting vector) over the waveguide's cross-section (e.g., plate
% thickness).
% 
% See also: powerFluxAxial, powerFluxTransverse
% 
% 2022-2024 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France

if ~isscalar(gew) % compute recursively for every waveguide problem in the vector "gew"
    compute = @(gewObj,datObj) powerFlux(gewObj, datObj); % function to apply
    P = arrayfun(compute,gew,dat,'UniformOutput',false); % apply to every object in the arrays "gew" and "dat"
    return;
end

p = poyntingVec(gew, dat);
P = GEWintegrate(gew, p);

end
