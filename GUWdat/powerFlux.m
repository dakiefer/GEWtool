function [P] = powerFlux(wguide, dat)
%POWERFLUX Compute total power flux along the waveguide.

p = poyntingVec(wguide, dat);

P = zeros(size(dat.k)); 
for i = 1:wguide.geom.nLay
    px = p{i}(:,:,:,1);             % power flux density along waveguide;
    lims = wguide.geom.yItf(i,:);   % lower and upper integration limits
    P = P + chebintegrate(px, lims, 3); % along 3rd dimension
end

end
