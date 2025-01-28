function [Px] = powerFluxAxial(dat)
% POWERFLUXAXIAL - Axial component of the power flux, i.e., along the wave vector.
%
% Computes the total power flux in direction of the wave vector (x-direction).
% This is done by integrating the power flux density's x-component (poynting
% vector component) over the waveguide's cross-section (thickness).
% 
% Usage: 
% > Px = powerFluxAxial(dat);
% 
% Arguments: 
% - dat:    GEWdat object that stores the computed modes. 
%
% See also: powerFlux, powerFluxTransverse
% 
% 2022-2024 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France

if ~isscalar(dat) % compute recursively for every waveguide problem in the vector "dat"
    Px = arrayfun(@powerFluxAxial,dat,'UniformOutput',false); % apply to every object in the arrays "dat"
    return;
end

if isa(dat.gew,"CylinderCircumferential")
    warning('Circumferential waves do not support this function yet. The results will be wrong.');
end
px = poyntingVecAxial(dat);
Px = GEWintegrate(dat.gew, px);

end
