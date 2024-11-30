function [py] = poyntingVecTransverse(dat)
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

if ~isscalar(dat) % compute recursively for every waveguide problem in the vector "dat"
    py = arrayfun(@poyntingVecTransverse,dat,'UniformOutput',false); % apply to every object in the arrays "dat"
    return;
end

v = velocity(dat);
T = stress(dat);     % inefficient: we compute components that we don't need
dof = 1:length(dat.gew.udof); % polarization

py = cell(dat.gew.geom.nLay,1);
for l = 1:dat.gew.geom.nLay
    if dat.gew.geom.Nudof(l) == 3
        py{l} = -1/2*sum(real(conj(v{l}).*T{l}(:,:,:,dof,2)), 4); % reduce T to components that yield py
    else
        py{l} = zeros([size(dat.w), dat.gew.geom.N(l)]);
    end
end

% TODO: at the moment, all layers are required to have the same number of
% polarization degrees of freedom, i.e., Nudof is the same for all layers.

end
