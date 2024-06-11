function [ceMag] = energyVel(gew, dat)
% energyVel - Alias to energyVelMag(). Energy velocity magnitude |ce|.
% 
% The energy velocity vector ce is the radio of the total power flux vector P to the 
% total stored energy H: ce = P/H . energyVel() is an alias to energyVelMag(), which
% returns the magnitude of ce.
% 
% See also: energyVelMag, energyVelVec
% 
% 2024 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France

ceMag = energyVelMag(gew, dat);

end
