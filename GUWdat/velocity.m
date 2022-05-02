function [v] = velocity(dat)
%VELOCITY compute the particle velocities

w = dat.w;
u = dat.u;
v = -1i*w.*squeeze(u);

end
