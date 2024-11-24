function [py] = poyntingVecTransverse(gew, dat)
% poyntingVecTransverse - Transverse component of the power flux density vectors. 
% 
% The Poynting vector is p = -1/2 v^* . T . Its real part represents the
% time-averaged power flux density vector. The transverse component is the
% y-component and it is orthogonal to the wave vector k (which is in ex direction) 
% and the plate's normal ez.
% 
% see also: poyntingVec, powerFlux
% 
% 2024 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France

if ~isscalar(gew) % compute recursively for every waveguide problem in the vector "gew"
    compute = @(gewObj,datObj) poyntingVecTransverse(gewObj, datObj); % function to apply
    py = arrayfun(compute,gew,dat,'UniformOutput',false); % apply to every object in the arrays "gew" and "dat"
    return;
end

v = velocity(dat);
T = stress(gew, dat);     % inefficient: we compute components that we don't need
dof = 1:length(gew.udof); % polarization

py = cell(gew.geom.nLay,1);
for l = 1:gew.geom.nLay
    if gew.geom.Nudof(l) == 3
        py{l} = -1/2*sum(real(conj(v{l}).*T{l}(:,:,:,dof,2)), 4); % reduce T to components that yield py
    else
        py{l} = zeros([size(dat.w), gew.geom.N(l)]);
    end
end

% TODO: at the moment, all layers are required to have the same number of
% polarization degrees of freedom, i.e., Nudof is the same for all layers.

end
