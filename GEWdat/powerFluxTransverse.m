function [Pz] = powerFluxTransverse(gew, dat)
% POWERFLUXTRANSVERSE - Compute the power flux transverse to the waveguide.
% Computes the power flux orthogonal to the wave vector (z-direction). This 
% is done by integrating the power flux density's z-component (poynting vector) 
% over the waveguide's cross-section (thickness). The transverse power flux 
% will usually be nonzero for propagation in nonprincipal directions of an 
% anisotropic material.
%
% See also: POWERFLUXAXIAL.
% 
% 2022-2023 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France

p = poyntingVec(gew, dat);

pz = cell(1, gew.geom.nLay); % allocate
for l = 1:gew.geom.nLay
    pz{l} = p{l}(:,:,:,3);        % power flux density along waveguide;
end

Pz = GEWintegrate(gew, pz);

end
