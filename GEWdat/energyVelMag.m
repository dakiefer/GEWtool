function [ceMag] = energyVelMag(dat)
% energyVelMag - Energy velocity magnitude |ce|.
% 
% The energy velocity vector ce is the ratio of the total power flux vector P to the 
% total stored energy H: ce = P/H . energyVelMag() returns the magnitude of ce.
% 
% See also: energyVel, energyVelAxial, energyVelTransverse, powerFlux
% 
% 2024 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France

if ~isscalar(dat) % compute recursively for every waveguide problem in the vector "dat"
    ceMag = arrayfun(@energyVelMag,dat,'UniformOutput',false); % apply to every object in the arrays "dat"
    return;
end

cex = energyVelAxial(dat);      % ex-direction
cey = energyVelTransverse(dat); % ey-direction
ceMag = sqrt(cex.^2 + cey.^2);

% NOTE: The above computation is done in terms of energy velocity components
% instead of power flux components on purpose because the former are better
% scaled (closer to unity) than the latter.

end
