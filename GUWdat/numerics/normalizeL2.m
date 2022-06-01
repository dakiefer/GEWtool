function uNorm = normalizeL2(u, w)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

uAbs = sum(conj(u).*u,4);
w = shiftdim(w(:), -2);

I = sum(w.*uAbs, 3);
uNorm = u./sqrt(I);

end