function uNorm = normalizeL2(u, w)
% normalizeL2 - normalize to L2-norm
% 
% Normalize the data u to unit L2-norm using the integration weights w. 
% The integration is performed along the 3rd dimension of u.
% Arguments:
% - u: data as array of dimension p >= 3
% - w: vector of integration weights
% Return values:
% - uNorm: array of dimension p-1
% 
% 2022 - Daniel A. Kiefer
% Institut Langevin, Paris, France

uAbs = sum(conj(u).*u,4);
w = shiftdim(w(:), -2);

I = sum(w.*uAbs, 3);
uNorm = u./sqrt(I);

end