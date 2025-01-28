function [Py] = powerFluxTransverse(dat)
% powerFluxTransverse - Transverse component of the power flux, orthogonal to k.
% 
% Computes the total power flux in z-direction. This component is orthogonal to
% the wave vector k (ex direction) and the plate's normal ey. This is done by
% integrating the power flux density's z-component (poynting vector component)
% over the waveguide's cross-section (thickness). The transverse power flux 
% will usually be nonzero for propagation along a non-symmetry direction of an 
% anisotropic material.
% 
% Usage: 
% > Py = powerFluxTransverse(dat);
% 
% Arguments: 
% - dat:    GEWdat object that stores the computed modes. 
%
% See also: powerFlux, powerFluxAxial
% 
% 2022-2024 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France

if ~isscalar(dat) % compute recursively for every waveguide problem in the vector "dat"
    Py = arrayfun(@powerFluxTransverse,dat,'UniformOutput',false); % apply to every object in the arrays "dat"
    return;
end

if isa(dat.gew,"CylinderCircumferential")
    warning('Circumferential waves do not support this function yet. The results will be wrong.');
end
py = poyntingVecTransverse(dat);
Py = GEWintegrate(dat.gew, py);

end
