function [pz] = poyntingVecTransverse(gew, dat)
% poyntingVecTransverse - Transverse component of the power flux density vectors. 
% 
% The Poynting vector is p = -1/2 v^* . T . Its real part represents the
% time-averaged power flux density vector. The transverse component is the
% z-component and it is orthogonal to the wave vector k (in ex direction) and
% the plate's normal ey.
% 
% see also: poyntingVec, powerFlux
% 
% 2024 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France

v = velocity(dat);
T = stress(gew, dat); % inefficient: we compute components that we don't need

pz = cell(gew.geom.nLay,1);
for l = 1:gew.geom.nLay
    pz{l} = -1/2*sum(real(conj(v{l}).*T{l}(:,:,:,:,3)), 4); % reduce T to components that yield pz
end

end
