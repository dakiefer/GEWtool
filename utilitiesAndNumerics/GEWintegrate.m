function [I] = GEWintegrate(gew, f, n)
% GEWintegrate - Integrate over all layers of the waveguide.
% 
% Integrates the data f along dimension n using the weights and limits 
% provided by gew. Arguments:
% - gew:  (Waveguide) Description of the waves
% - f:    (cell array) Data to be integrated. f{1} corresponds to the data of
%         the first layer and is a p-dimensional array.
% - n:    (integer, default: 3) Dimension along which to integrate. The third
%         dimension of the array f{1} will usually correspond to samples on the
%         waveguide cross-section.
% 
% Return value: 
% - I:    (array of dimension p-1) Integrated values.
%
% 2022-2024 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France

if nargin < 3 
    n = 3; % dimension to be integrated;
end
s = size(f{1}); s(n) = []; % remove dimension to be integrated
I = zeros(s);
for l = 1:gew.geom.nLay  % for every layer
    lay = gew.lay(l);
    hi = gew.geom.hl(l); % thickness of layer i
    w = shiftdim(lay.w(:), -n+1);       % integration weights moved to dim n
    if isa(gew,"Cylinder") && lay.r(1) ~= 0
        r = shiftdim(lay.r(:), -n+1); % shifted to dimension of integration
        integrand = f{l}.*r;
    elseif isa(gew,"Cylinder") && lay.r(1) == 0 % avoid integrating at r = 0
        [flNew, ri, wi] = interpolateToGaussLegendre(lay, f{l}, n);
        w = wi; % replace integration weights with those of the new grid
        integrand = flNew.*ri;
    else
        integrand = f{l};
    end
    I = I + reshape( sum(w.*integrand*hi, n) , s ); % reshape removes the singleton dimension
end

end

function [flNew, ri, wi] = interpolateToGaussLegendre(lay, fl, n)
    dimsFl = ndims(fl);
    [yi, wi] = lgwt(lay.N,0,1); % Gauss-Legendre points and weights
    yi = flip(yi); % upside-down
    wi = shiftdim(wi, -n+1);      % replace by integration weights on integration nodes yi
    % interpolate each component of f onto new coordinates yi:
    fl = shiftdim(fl,(n-1)); % integration on first dimension
    fNorm = norm(fl,'fro');
    comps = reshape(fl, lay.N, [])/fNorm; % all components of the data
    for k = 1:size(comps,2)
        data = [lay.eta(:), comps(:,k)];
        comps(:,k) = barylag(data,yi); % replace with interpolated data
    end
    flNew = shiftdim(reshape(comps, size(fl)), dimsFl-n+1)*fNorm; % positive shift is performed circularly while negative ones add singleton dimensions
    ri = shiftdim(yi*lay.h, -n+1);    % scale to SI units, inner radius is zero
end