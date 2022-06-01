function [I] = GEWintegrate(wguide, f, n)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

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
