function [px] = poyntingVecAxial(gew, dat)
% poyntingVecAxial - Axial component of the power flux density vectors. 
% 
% The Poynting vector is p = -1/2 v^* . T . Its real part represents the
% time-averaged power flux density vector. The x-component is aligned with the
% wave vector k and is the "axial component". 
% 
% see also: poyntingVec, velocity, stress, powerFlux
% 
% 2024 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France

if ~isscalar(gew) % compute recursively for every waveguide problem in the vector "gew"
    compute = @(gewObj,datObj) poyntingVecAxial(gewObj, datObj); % function to apply
    px = arrayfun(compute,gew,dat,'UniformOutput',false); % apply to every object in the arrays "gew" and "dat"
    return;
end

v = velocity(dat);
T = stress(gew, dat); % inefficient: we compute components that we don't need
udof = gew.udof;      % polarization

px = cell(gew.geom.nLay,1);
for l = 1:gew.geom.nLay
    px{l} = -1/2*sum(real(conj(v{l}).*T{l}(:,:,:,udof,1)), 4); % reduce T to components that yield px
end

end
