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
dof = 1:length(dat.gew.udof);      % polarization

p = cell(dat.gew.geom.nLay,1);
for l = 1:dat.gew.geom.nLay
    pl = -1/2*sum(real(conj(v{l}).*T{l}(:,:,:,dof,:)), 4);
    p{l} = permute(pl, [1 2 3 5 4]);
end

end
