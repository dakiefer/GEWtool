function [ce] = energyVelVec(gew, dat)
% ENERGYVEL - Energy velocity vector ce.
% The energy velocity vector ce is the radio of the total power flux vector P to the 
% total stored energy H: ce = P/H . 
% 
% See also: energyVelAxial, energyVelTransverse, powerFlux
% 
% 2024 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France

P = powerFlux(gew, dat);   % for a plate: vector in the plane of the plate
H = energyTotal(gew, dat);
ce = P./H;

% NOTE: this function depends on many of the field calculation functions:
% - GEWintegrate
% - powerFlux
%   - poyntingVec
%       - velocity
%       - stress
% - energyTotal 
%   - energyKinetic
%       - energyDensityKinetic
%           - velocity
% 
% Additionally, without equipartion of energy, it also epends on
%   - energyElastic
%       - energyDensityElastic
%           - strain
%           - stress
% 
% Accordingly, energy velocity computations can be used to test the correct
% impementation of the field computation functions.

end
