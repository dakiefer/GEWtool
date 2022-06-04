function [P] = powerFlux(wguide, dat)
%POWERFLUX Compute total power flux along the waveguide.

p = poyntingVec(wguide, dat);
for i = 1:wguide.geom.nLay
    px{i} = p{i}(:,:,:,1);        % power flux density along waveguide;
end

P = GEWintegrate(wguide, px);

end
