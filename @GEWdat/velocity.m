function [v] = velocity(dat)
%VELOCITY compute the particle velocities
% 
% Usage: 
% > v = velocity(dat);
% 
% Arguments: 
% - dat:    GEWdat object that stores the computed modes. 
%

if ~isscalar(dat) % compute recursively for every waveguide problem in the vector "gew"
    v = arrayfun(@velocity,dat,'UniformOutput',false); % apply to every object in the array "dat"
    return;
end

w = dat.w;
u = dat.u;
v = cell(size(u));
for i = 1:length(v) % for every layer
    v{i} = -1i*w.*u{i};
end

end
