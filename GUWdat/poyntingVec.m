function [p] = poyntingVec(wguide, dat)
%POYNTINGVEC compute the poynting vectors, i.e., power flux density vectors. 

v = velocity(dat);
T = stress(wguide, dat);
p = -1/2*real(squeeze(sum(conj(v).*T, 4)));

end