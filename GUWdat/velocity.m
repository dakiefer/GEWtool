function [v] = velocity(dat)
%VELOCITY compute the particle velocities

w = dat.w;
u = dat.u;
v = cell(size(u));
for i = 1:length(v) % for every layer
    v{i} = -1i*w.*u{i};
end

end
