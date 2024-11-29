function [T] = stress(dat)
% stress - Stress tensor T.
% 
% Computed by double-contraction of the displacement gradient F with the
% stiffness tensor c, i.e.,
% 
% T = c : F
% 
% 2024 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France

if ~isscalar(dat) % compute recursively for every waveguide problem in the vector "dat"
    T = arrayfun(@stress,dat,'UniformOutput',false); % apply to every object in the arrays "dat"
    return;
end

if isa(dat.gew, 'Plate')
    udof = dat.gew.udof; % polarization (which of the ux, uy, uz components are represented)
else
    udof = 1:3; % Cylinder always has all stresses
end
F = displGrad(dat.gew, dat);
T = cell(dat.gew.geom.nLay, 1); % allocate for each layer
for l = 1:dat.gew.geom.nLay
    c = dat.gew.lay{l}.mat.c; % stiffness tensor
    c = shiftdim(c(udof,udof,udof,udof), -3); % shift to the correct dimension
    Fl = permute(F{l}, [1 2 3 6 7 4 5]); % additional dimension for contraction with c
    T{l} = sum(sum(c.*Fl, 6), 7); % double contraction 
end

end
