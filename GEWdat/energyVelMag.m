function [ceMag] = energyVelMag(gew, dat)
% energyVelMag - Energy velocity magnitude |ce|.
% 
% The energy velocity vector ce is the radio of the total power flux vector P to the 
% total stored energy H: ce = P/H . energyVelMag() returns the magnitude of ce.
% 
% See also: energyVel, energyVelAxial, energyVelTransverse, powerFlux
% 
% 2024 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France

cex = energyVelAxial(gew, dat);      % ex-direction
cez = energyVelTransverse(gew, dat); % ez-direction
ceMag = sqrt(cex.^2 + cez.^2);

% NOTE: The above computation is done in terms of energy velocity components
% instead of power flux components on purpose because the former are better
% scaled (closer to unity) than the latter.

end