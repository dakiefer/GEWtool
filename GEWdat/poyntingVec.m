function [p] = poyntingVec(gew, dat)
% poyntingVec - Power flux density vectors.
% 
% The Poynting vector is p = -1/2 v^* . T . Its real part represents the
% time-averaged power flux density vectors. 
% 
% See also: poyntingVec, velocity, stress, powerFlux
% 
% 2024 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France

v = velocity(dat);
T = stress(gew, dat);

p = cell(gew.geom.nLay,1);
for i = 1:gew.geom.nLay
    p_i = -1/2*sum(real(conj(v{i}).*T{i}), 4);
    p{i} = permute(p_i, [1 2 3 5 4]);
end

end
