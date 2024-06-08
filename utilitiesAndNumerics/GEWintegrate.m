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

if isa(gew,"Cylinder")
    warning('Cylinders do not support this function yet. The results will be wrong.');
end
if nargin < 3 
    n = 3; % dimension to be integrated;
end
s = size(f{1}); s(n) = []; % remove dimension to be integrated
I = zeros(s);
for i = 1:gew.geom.nLay  % for every layer
    hi = gew.geom.hl(i); % thickness of layer i
    w = shiftdim(gew.lay(i).w(:), -n+1);       % integration weights moved to dim n
    I = I + reshape( sum(w.*f{i}*hi, n) , s ); % reshape removes the singleton dimension
end

end
