function [T] = stress(gew, dat)
% stress - Stress tensor T.
% 
% Computed by double-contraction of the displacement gradient F with the
% stiffness tensor c, i.e.,
% 
% T = c : F
% 
% 2024 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France

if ~isscalar(gew) % compute recursively for every waveguide problem in the vector "gew"
    compute = @(gewObj,datObj) stress(gewObj, datObj); % function to apply
    T = arrayfun(compute,gew,dat,'UniformOutput',false); % apply to every object in the arrays "gew" and "dat"
    return;
end

if isa(gew, 'Plate')
    udof = gew.udof; % polarization (which of the ux, uy, uz components are represented)
else
    udof = 1:3; % Cylinder always has all stresses
end
F = displGrad(gew, dat);
T = cell(gew.geom.nLay, 1); % allocate for each layer
for l = 1:gew.geom.nLay
    c = gew.lay(l).mat.c; % stiffness tensor
    c = shiftdim(c(udof,udof,udof,udof), -3); % shift to the correct dimension
    Fl = permute(F{l}, [1 2 3 6 7 4 5]); % additional dimension for contraction with c
    T{l} = sum(sum(c.*Fl, 6), 7); % double contraction 
end

end
