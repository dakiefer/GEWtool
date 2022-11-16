function [P] = powerFlux(gew, dat)
%POWERFLUX Compute total power flux along the waveguide.

p = poyntingVec(gew, dat);

px = cell(1, gew.geom.nLay); % allocate
for l = 1:gew.geom.nLay
    px{l} = p{l}(:,:,:,1);        % power flux density along waveguide;
end

P = GEWintegrate(gew, px);

end
