function [T] = stress(dat)
% stress - Stress tensor T.
% 
% Mechanics: computed by double-contraction of the displacement gradient F with the
% stiffness tensor c, i.e.,
% 
% T = c : F 
% 
% Piezoelectrics: contraction of the displacement gradient F and potential gradient G
% with the stiffness tensor c and piezoelectric stress constants e, i.e.,
% 
% T = c : F + G . e
% 
% Usage: 
% > T = stress(dat);
% 
% Arguments: 
% - dat:    GEWdat object that stores the computed modes. 
%
% 2024-2026 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France

if ~isscalar(dat) % compute recursively for every waveguide problem in the vector "dat"
    T = arrayfun(@stress,dat,'UniformOutput',false); % apply to every object in the array "dat"
    return;
end

if isa(dat.gew, 'Plate')
    udof = dat.gew.udof; % polarization (which of the ux, uy, uz components are represented)
else
    udof = 1:3; % Cylinder always has all stresses
end

% retrieve required field(s): 
F = displGrad(dat.gew, dat);
if isPiezoelectric(dat.gew)
    G = potentialGrad(dat.gew, dat); 
end

T = cell(dat.gew.geom.nLay, 1); % allocate for each layer
for l = 1:dat.gew.geom.nLay
    matl = dat.gew.lay{l}.mat; 
    if isa(matl, "MaterialPiezoelectric")
        T{l} = matl.stress(F{l}, G{l}, dat.gew.np, udof);
    else
        T{l} = matl.stress(F{l}, dat.gew.np, udof);
    end
end

% for l = 1:dat.gew.geom.nLay
%     cl = shiftdim( dat.gew.lay{l}.mat.c(udof,udof,udof,udof), -3)/dat.gew.np.c0; % stiffness: shift to correct dimension
%     el = shiftdim( dat.gew.lay{l}.mat.e(udof,udof,udof), -3)/dat.gew.np.e0; % piezoelec. coupling: shift to correct dimension
%     Fl = permute(F{l}, [1 2 3 6 7 4 5]); % additional dimensions 6 and 7 for contraction with c
%     El = E{l};
%     T{l} = sum(sum(cl.*Fl, 6), 7) - permute( sum( El.*el , 4), [1 2 3 5 6 4] ); % permute -> remove contracted dimension
% end

end
