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

p = poyntingVec(gew, dat);
P = GEWintegrate(gew, p);

end
