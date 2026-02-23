function [p] = poyntingVec(dat)
% poyntingVec - Power flux density vectors.
% 
% The Poynting vector is p = -1/2 v^* . T . Its real part represents the
% time-averaged power flux density vectors. 
% 
% Usage: 
% > p = poyntingVec(dat);
% 
% Arguments: 
% - dat:    GEWdat object that stores the computed modes. 
%
% See also: poyntingVec, velocity, stress, powerFlux
% 
% 2024 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France

if ~isscalar(dat) % compute recursively for every waveguide problem in the vector "dat"
    p = arrayfun(@poyntingVec,dat,'UniformOutput',false); % apply to every object in the array "dat"
    return;
end

v = velocity(dat);
T = stress(dat);
udof = dat.gew.udof;
if size(v{1},4) == size(T{1},4) % reduced plain strain, i.e., [ux, uz] -> indices [1 2] instead of [1 3]
    udof = 1:length(udof);
end

if isPiezoelectric(dat.gew)
    phi = potential(dat); 
    D   = electricFluxDensity(dat); 
    w   = dat.w/dat.gew.np.fh0*dat.gew.np.h0; % normalized as fields (phi, v, etc) 
end

p = cell(dat.gew.geom.nLay,1);
for l = 1:dat.gew.geom.nLay
    p{l} = -1/2*sum(real( conj(v{l}(:,:,:,udof)).*T{l}(:,:,:,udof,:) ), 4); % except for dof, v = 0
    p{l} = permute(p{l}, [1 2 3 5 4]); % removes 4th dimension (singleton) 
    if isPiezoelectric(dat.gew)
        plElec =  1/2*real( phi{l}.*conj(-1i*w.*D{l}(:,:,:,udof)) ) ;
        p{l}(:,:,:,udof) = p{l}(:,:,:,udof) + plElec(:,:,:,udof); % plElec might be 2x1 while plMech could be 3x1 
    end
end

end
