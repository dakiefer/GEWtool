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

if ~isscalar(gew) % compute recursively for every waveguide problem in the vector "gew"
    compute = @(gewObj,datObj) poyntingVecTransverse(gewObj, datObj); % function to apply
    pz = arrayfun(compute,gew,dat,'UniformOutput',false); % apply to every object in the arrays "gew" and "dat"
    return;
end

v = velocity(dat);
T = stress(gew, dat); % inefficient: we compute components that we don't need
udof = gew.udof;      % polarization

pz = cell(gew.geom.nLay,1);
for l = 1:gew.geom.nLay
    if gew.geom.Nudof(l) == 3
        pz{l} = -1/2*sum(real(conj(v{l}).*T{l}(:,:,:,udof,3)), 4); % reduce T to components that yield pz
    else
        pz{l} = zeros([size(dat.w), gew.geom.N(l)]);
    end
end

% TODO: at the moment, all layers are required to have the same number of
% polarization degrees of freedom, i.e., Nudof is the same for all layers.

end
