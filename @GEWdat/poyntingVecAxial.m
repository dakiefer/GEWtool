function [px] = poyntingVecAxial(dat)
% poyntingVecAxial - Axial component of the power flux density vectors. 
% 
% The Poynting vector is p = -1/2 v^* . T . Its real part represents the
% time-averaged power flux density vector. The x-component is aligned with the
% wave vector k and is the "axial component". 
% 
% see also: poyntingVec, velocity, stress, powerFlux
% 
% 2024 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France

if ~isscalar(dat) % compute recursively for every waveguide problem in the vector "dat"
    px = arrayfun(@poyntingVecAxial,dat,'UniformOutput',false); % apply to every object in the arrays "dat"
    return;
end

v = velocity(dat);
T = stress(dat); % inefficient: we compute components that we don't need
x = dat.gew.udofAxial; % 1 for Plate and Cylinder, 2 for Circumferential

px = cell(dat.gew.geom.nLay,1);
for l = 1:dat.gew.geom.nLay
    if size(v{l},4) == size(T{l},4)
        px{l} = -1/2*sum(real(conj(v{l}).*T{l}(:,:,:,:,x)), 4); % reduce T to components that yield px
    else
        dof = dat.gew.udof; 
        px{l} = -1/2*sum(real(conj(v{l}).*T{l}(:,:,:,dof,x)), 4); % except for dof, v = 0
    end
end

end
