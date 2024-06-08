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

px = poyntingVecAxial(gew,dat);
Px = GEWintegrate(gew, px);

end
