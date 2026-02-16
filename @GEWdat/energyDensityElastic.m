function [eElastic] = energyDensityElastic(dat)
%ENERGYDENSITYELASTIC Compute the elastic energy density.
% 
% eElastic = 1/4*real(S':c:S)
%
% Usage: 
% > eElastic = energyDensityElastic(dat);
% 
% Arguments: 
% - dat:    GEWdat object that stores the computed modes. 
%

if ~isscalar(dat) % compute recursively for every waveguide problem in the vector "dat"
    eElastic = arrayfun(@energyDensityElastic,dat,'UniformOutput',false); % apply to every object in the array "dat"
    return;
end

if isDissipative(dat.gew)
    warning('energyDensityElastic(): Is it valid for viscous media?')
end

S = strain(dat);
if isa(dat.gew, 'Plate')
    udof = dat.gew.udof; % polarization (which of the ux, uy, uz components are represented)
else
    udof = 1:3; % Cylinder always has all stresses
end

eElastic = cell(dat.gew.geom.nLay,1);
for l = 1:dat.gew.geom.nLay
    c = dat.gew.lay{l}.mat.c/dat.gew.np.c0; % stiffness tensor
    c = shiftdim(c(udof,udof,udof,udof), -3); % shift to the correct dimension
    Sl = permute(S{l}, [1 2 3 6 7 4 5]); % additional dimension for contraction with c
    Tl = sum(sum(c.*Sl, 6), 7); % avoid colling stress(): T is not the same as it includes piezoelectric coupling
    eElastic{l} = 1/4*real(sum(sum(conj(S{l}).*Tl, 5), 4)); % double contraction
end

end
