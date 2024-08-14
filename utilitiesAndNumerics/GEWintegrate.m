function [I] = GEWintegrate(gew, f, ls, dim)
% GEWintegrate - Integrate over all (or specified) layers of the waveguide.
% 
% Integrates the data "f" along dimension "dim" using the quadrature nodes and
% weights provided by the Waveguide object "gew". The limits of integration are
% also determined from "gew". For a cylindrical layer beginning at r = 0 the
% field is first interpolated to Gauss-Legendre points to avoid the singularity.
% 
% Arguments:
% - gew:  (Waveguide) Description of the waves
% - f:    (cell array) Data to be integrated. f{1} corresponds to the data of
%         the first layer and is a p-dimensional array.
% - ls:   (integer vector, default: 1:gew.geom.nLay) layer numbers to include in 
%         integration. You can specify just one layer number or several. If
%         empty, all layers are integrated per default. 
% - dim:  (integer, default: 3) Dimension of the numeric array along which to 
%         integrate. The third dimension of the array f{1} will usually correspond 
%         to samples on the waveguide cross-section.
% 
% Return value: 
% - I:    (array of dimension p-1) Integrated values.
% 
% See also: GEWintegrateEachLayer, lglnodes, lgwt
%
% 2022-2024 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France

if nargin < 4
    dim = 3; % dimension to be integrated;
end
if nargin < 3 || isempty(ls)
    ls =  1:gew.geom.nLay;
end
s = size(f{1}); s(dim) = []; % remove dimension to be integrated
if length(s) < 2, s(end+1) = 1; end % matlab needs strange 2d-size for vectors
I = zeros(s);
for l = ls  % for every layer specified in ls
    lay = gew.lay{l};
    hi = gew.geom.hl(l); % thickness of layer i
    w = shiftdim(lay.w(:), -dim+1);       % integration weights moved to dim n
    if isa(gew,"Cylinder") && lay.r(1) ~= 0
        r = shiftdim(lay.r(:), -dim+1); % shifted to dimension of integration
        integrand = f{l}.*r;
    elseif isa(gew,"Cylinder") && lay.r(1) == 0 % avoid integrating at r = 0
        [flNew, ri, wi] = interpolateToGaussLegendre(lay, f{l}, dim);
        w = wi; % replace integration weights with those of the new grid
        integrand = flNew.*ri;
    else
        integrand = f{l};
    end
    I = I + reshape( sum(w.*integrand*hi, dim) , s ); % reshape removes the singleton dimension
end

end

function [flNew, ri, wi] = interpolateToGaussLegendre(lay, fl, n)
    dimsFl = ndims(fl);
    [yi, wi] = lgwt(lay.N,0,1); % Gauss-Legendre points and weights
    yi = flip(yi); % upside-down
    wi = shiftdim(wi, -n+1);      % replace by integration weights on integration nodes yi
    % interpolate each component of f onto new coordinates yi:
    fl = shiftdim(fl,(n-1)); % integration on first dimension
    fNorm = norm(fl(:));
    comps = reshape(fl, lay.N, [])/fNorm; % all components of the data
    for k = 1:size(comps,2)
        data = [lay.eta(:), comps(:,k)];
        comps(:,k) = barylag(data,yi); % replace with interpolated data
    end
    flNew = shiftdim(reshape(comps, size(fl)), dimsFl-n+1)*fNorm; % positive shift is performed circularly while negative ones add singleton dimensions
    ri = shiftdim(yi*lay.h, -n+1);    % scale to SI units, inner radius is zero
end