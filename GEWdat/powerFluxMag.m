function [Pmag] = powerFluxMag(gew, dat)
% powerFluxMag - Power flux magnitude.
%
% Computes the magnitude of the power flux vectors P. In a plate, the power flux
% vector P might ly in any direction within the plane of the plate.
% 
% See also: powerFlux, powerFluxAxial, powerFluxTransverse
% 
% 2022-2024 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France

if ~isscalar(gew) % compute recursively for every waveguide problem in the vector "gew"
    compute = @(gewObj,datObj) powerFluxMag(gewObj, datObj); % function to apply
    Pmag = arrayfun(compute,gew,dat,'UniformOutput',false); % apply to every object in the arrays "gew" and "dat"
    return;
end

normP = gew.np.c0 * gew.np.fh0 * gew.np.h0;
Px = powerFluxAxial(gew, dat)/normP;
Pz = powerFluxTransverse(gew, dat)/normP; 
Pmag  = sqrt(Px.^2 + Pz.^2)*normP; 

% NOTE: The power flux in SI units has very high values. Do the computation in
% normalized values and scale back. 

end
