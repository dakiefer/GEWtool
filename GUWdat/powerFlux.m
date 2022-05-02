function [P] = powerFlux(wguide, dat)
%POWERFLUX Compute total power flux along the waveguide.

p = poyntingVec(wguide, dat);
px = p(:,:,:,1); % along waveguide;
P = chebintegrate(px, wguide.geom.yItf, 3);

end
