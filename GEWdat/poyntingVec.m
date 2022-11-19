function [p] = poyntingVec(gew, dat)
%POYNTINGVEC compute the poynting vectors, i.e., power flux density vectors. 

v = velocity(dat);
T = stress(gew, dat);

p = cell(gew.geom.nLay,1);
for i = 1:gew.geom.nLay
    p_i = -1/2*sum(real(conj(v{i}).*T{i}), 4);
    p{i} = permute(p_i, [1 2 3 5 4]);
end

end
