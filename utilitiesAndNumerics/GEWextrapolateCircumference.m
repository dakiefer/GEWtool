function [ue, phi] = GEWextrapolateCircumference(u, n, Ntheta)
% GEWextrapolateCircumference - Extrapolate data along the angular direction.
% 
% The modal field of waves in a cylinder are 
%
%    u(x,φ,r) = u(r)*exp(i*k*x)*exp(i*n*φ) 
% 
% for all (displacement or other) components of the vector/tensor u. Note that
% the tensor u(r) implicitly depends on φ through the cylindrical basis vectors
% and this is taken into account by GEWtool during the computation. However, the
% "j"th component uj(r) is independent of φ.
% 
% The uj(r) have been computed by computeK() or computeW(). GEWextrapolateCircumference
% performs the extrapolation with exp(i*n*φ), where "n" is the circumferential
% wavenumber and needs to be provided as an argument.
% 
% Arguments: 
% - u:      (array of dim >= 1) field data of _one_ selected mode to be 
%           extrapolated, e.g., ''dat.u{1}(1,1,:,:)'' returned by computeW() or 
%           computeK(). GEWinterpolate() directly generates data 'u' in the
%           appropriate structure.
% - n:      (scalar) circumferential order (as prescribed when choosing the wave family)
% - Ntheta: (scalar, default: 60) number of points along the circumference.
% 
% Return values: 
% - ue:     (array) extrapolated field data. It has one more dimension than 'u' that
%           describes the angular dependence. The new dimension is the second one.
%           dimension (the first one is usually the radial dependance).
% - theta:  (vector) vector of angular coordinates on which the extrapolation was
%           performed.
%
% See also: GEWinterpolate, extractModes.
% 
% 2024 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France

if nargin < 3
    Ntheta = 60;    % default 
end
if iscell(u)        % does not support cell array
    error('GEWTOOL:GEWextrapolateCircumference', ...
        'A cell array has been provided as field data. Extract the data into an array. You can do so by providing the data of only one layer or by interpolating onto a grid using GEWinterpolate().');
end
u = squeeze(u);     % ignore singleton dimensions
if isvector(u)
    nd = 1;         % one dimension
else
    nd = ndims(u);  % 'nd' number of dimensions
end
perm = [1, nd+1, 2:nd];          % to create a new dimension at second position
u = permute(u,perm);
phi = linspace(0, 2*pi, Ntheta); % angular coordinates
ue = u.*exp(1i*n*phi);           % extrapolate using circumf. order 'n'

end
