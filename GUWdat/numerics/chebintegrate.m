function [I] = chebintegrate(f, lim, dim)
% CHEBINTEGRATE Integrates the data f sampled on Chebyshev points
% from lim(1) to lim(2) along dimension dim.
% - f: array of function samples to be integrated
% - lim: 2-vector with lower and upper limit (in that order)
% - dim: dimension to integrate along (scalar, default: 1)

if nargin<3, dim=1; end
[~, w] = fclencurt(size(f,dim), lim(1), lim(2));
w = shiftdim(w, -dim+1);
I = sum(w.*f, dim);

end
