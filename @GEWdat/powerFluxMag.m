function [Pmag] = powerFluxMag(dat)
% powerFluxMag - Power flux magnitude.
%
% Computes the magnitude of the power flux vectors P. In a plate, the power flux
% vector P might ly in any direction within the plane of the plate.
% 
% See also: powerFlux, powerFluxAxial, powerFluxTransverse
% 
% 2022-2024 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France

if ~isscalar(dat) % compute recursively for every waveguide problem in the vector "gew"
    Pmag = arrayfun(@powerFluxMag,dat,'UniformOutput',false); % apply to every object in the arrays "dat"
    return;
end

if isa(dat.gew,"CylinderCircumferential")
    warning('Circumferential waves do not support this function yet. The results will be wrong.');
end
normP = dat.gew.np.c0 * dat.gew.np.fh0 * dat.gew.np.h0;
Px = powerFluxAxial(dat)/normP;
Py = powerFluxTransverse(dat)/normP; 
Pmag  = sqrt(Px.^2 + Py.^2)*normP; 

% NOTE: The power flux in SI units has very high values. Do the computation in
% normalized values and scale back. 

end
