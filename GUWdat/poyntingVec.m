function [p] = poyntingVec(wguide, dat)
%POYNTINGVEC compute the poynting vectors, i.e., power flux density vectors. 

v = velocity(dat);
T = stress(wguide, dat);

p = cell(wguide.geom.nLay,1);
for i = 1:wguide.geom.nLay
    p{i} = -1/2*real(sum(conj(v{i}).*T{i}, 5));
end

end
