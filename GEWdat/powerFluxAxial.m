function [Px] = powerFluxAxial(gew, dat)
% POWERFLUXAXIAL - Compute the total power flux along the waveguide.
% Computes the power flux in direction of the wave vector (x-direction). This 
% is done by integrating the power flux density's x-component (poynting vector) 
% over the waveguide's cross-section (thickness).
% 
% See also: POWERFLUXTRANSVERSE.
% 
% 2022-2023 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France

px = poyntingVecAxial(gew,dat);
Px = GEWintegrate(gew, px);

end
