function [p] = poyntingVec(dat)
% poyntingVec - Power flux density vectors.
% 
% The Poynting vector is p = -1/2 v^* . T . Its real part represents the
% time-averaged power flux density vectors. 
% 
% See also: poyntingVec, velocity, stress, powerFlux
% 
% 2024 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France

if ~isscalar(dat) % compute recursively for every waveguide problem in the vector "dat"
    p = arrayfun(@poyntingVec,dat,'UniformOutput',false); % apply to every object in the arrays "dat"
    return;
end

v = velocity(dat);
T = stress(dat);

p = cell(dat.gew.geom.nLay,1);
for l = 1:dat.gew.geom.nLay
    if size(v{l},4) == size(T{l},4)
        pl = -1/2*sum(real(conj(v{l}).*T{l}), 4);
    else
        dof = dat.gew.udof; 
        pl = -1/2*sum(real(conj(v{l}).*T{l}(:,:,:,dof,:)), 4); % except for dof, v = 0
    end
    p{l} = permute(pl, [1 2 3 5 4]); % removes 4th dimension (singleton) 
end

end
