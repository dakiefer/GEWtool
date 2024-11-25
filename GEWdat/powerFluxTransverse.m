function [Pz] = powerFluxTransverse(gew, dat)
% powerFluxTransverse - Transverse component of the power flux, orthogonal to k.
% 
% Computes the total power flux in z-direction. This component is orthogonal to
% the wave vector k (ex direction) and the plate's normal ey. This is done by
% integrating the power flux density's z-component (poynting vector component)
% over the waveguide's cross-section (thickness). The transverse power flux 
% will usually be nonzero for propagation along a non-symmetry direction of an 
% anisotropic material.
% 
% See also: powerFlux, powerFluxAxial
% 
% 2022-2024 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France

if ~isscalar(gew) % compute recursively for every waveguide problem in the vector "gew"
    compute = @(gewObj,datObj) powerFluxTransverse(gewObj, datObj); % function to apply
    Pz = arrayfun(compute,gew,dat,'UniformOutput',false); % apply to every object in the arrays "gew" and "dat"
    return;
end

if isa(gew,"CylinderCircumferential")
    warning('Circumferential waves do not support this function yet. The results will be wrong.');
end
pz = poyntingVecTransverse(gew, dat);
Pz = GEWintegrate(gew, pz);

end
