function [Px] = powerFluxAxial(gew, dat)
% POWERFLUXAXIAL - Axial component of the power flux, i.e., along the wave vector.
%
% Computes the total power flux in direction of the wave vector (x-direction).
% This is done by integrating the power flux density's x-component (poynting
% vector component) over the waveguide's cross-section (thickness).
% 
% See also: powerFlux, powerFluxTransverse
% 
% 2022-2024 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France

if ~isscalar(gew) % compute recursively for every waveguide problem in the vector "gew"
    compute = @(gewObj,datObj) powerFluxAxial(gewObj, datObj); % function to apply
    Px = arrayfun(compute,gew,dat,'UniformOutput',false); % apply to every object in the arrays "gew" and "dat"
    return;
end

px = poyntingVecAxial(gew,dat);
Px = GEWintegrate(gew, px);

end
