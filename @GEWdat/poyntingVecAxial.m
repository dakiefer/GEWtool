function [px] = poyntingVecAxial(dat)
% poyntingVecAxial - Axial component of the power flux density vectors. 
% 
% The Poynting vector is p = -1/2 v^* . T . Its real part represents the
% time-averaged power flux density vector. The x-component is aligned with the
% wave vector k and is the "axial component". 
% 
% Usage: 
% > px = poyntingVecAxial(dat);
% 
% Arguments: 
% - dat:    GEWdat object that stores the computed modes. 
%
% see also: poyntingVec, velocity, stress, powerFlux
% 
% 2024 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France

if ~isscalar(dat) % compute recursively for every waveguide problem in the vector "dat"
    px = arrayfun(@poyntingVecAxial,dat,'UniformOutput',false); % apply to every object in the array "dat"
    return;
end

v = velocity(dat);
T = stress(dat); % inefficient: we compute components that we don't need
x = dat.gew.udofAxial; % 1 for Plate and Cylinder, 2 for Circumferential
udof = dat.gew.udof;
if size(v{1},4) == size(T{1},4) % reduced plain strain, i.e., [ux, uz] -> indices [1 2] instead of [1 3]
    udof = 1:length(udof);
end

if isPiezoelectric(dat.gew)
    phi = potential(dat); 
    D   = electricFluxDensity(dat); 
    w   = dat.w/dat.gew.np.fh0*dat.gew.np.h0; % normalized as fields (phi, v, etc) 
end

px = cell(dat.gew.geom.nLay,1);
for l = 1:dat.gew.geom.nLay
    px{l} = -1/2*sum(real( conj(v{l}).*T{l}(:,:,:,udof,x) ), 4); % except for dof, v = 0
    if isPiezoelectric(dat.gew)
        plElec =  1/2*real( phi{l}.*conj(-1i*w.*D{l}(:,:,:,x)) ) ;
        px{l} = px{l} + plElec;
    end
end

end
