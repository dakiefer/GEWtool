function [I] = GEWintegrate(wguide, f, n)
% GEWintegrate - Integrate over all layers of the waveguide.
% 
% Integrates the data f along dimension n using the weights and limits 
% provided by wguide. Arguments:
% - wguide: description of the waveguide (class Waveguide) 
% - f: data to be integrated (array of dimension p >= 3)
% - n: dimension along which to integrate (integer, default: 3)
% Return value: 
% - I: Integrated values (array of dimension p-1)
%
% 2022 - Daniel A. Kiefer
% Institut Langevin, Paris, France

if nargin < 3 
    n = 3; % dimension to be integrated;
end
s = size(f{1}); s(n) = []; % remove dimension to be integrated
I = zeros(s);
for i = 1:wguide.geom.nLay % for every layer
    hi = wguide.geom.h(i); % thickness of layer i
    w = shiftdim(wguide.lay(i).w(:), -n+1); % integration weights moved to dim. n
    I = I + sum( w.*f{i}*hi , n);
end

end