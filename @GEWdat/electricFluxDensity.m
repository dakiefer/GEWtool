function [D] = electricFluxDensity(dat)
% electricFluxDensity - returns the electric flux density D. 
% 
% D = e:F - eps.G = e:grad(u) - eps.grad(phi)
% 
% Only relevant for piezoelectric waveguides. A zero-cell array of appropriate size 
% is returned otherwise. 
% 
% Usage: 
% > D = electricFluxDensity(dat)
% 
% Arguments: 
% - dat:    GEWdat object that stores the computed modes. 
%

if ~isscalar(dat) % compute recursively for every waveguide problem in the vector "gew"
    D = arrayfun(@electricFluxDensity,dat,'UniformOutput',false); % apply to every object in the array "dat"
    return;
end

G = potentialGrad(dat.gew,dat); % = zero if the plate is not piezoelectric
F = displGrad(dat.gew,dat); 
udof = dat.gew.udof; % E has same components as displacements u

D = cell(dat.gew.geom.nLay, 1); % allocate for each layer
for l = 1:dat.gew.geom.nLay
    matl = dat.gew.lay{l}.mat;
    if isa(matl, "MaterialPiezoelectric")
        D{l} = matl.electricFluxDensity(F{l}, G{l}, dat.gew.np, udof);
    else
        D{l} = G{l}; % = 0
    end
end

% for l = 1:dat.gew.geom.nLay
%     epsl = shiftdim( dat.gew.lay{l}.mat.epsilon(udof,udof), -3)/dat.gew.np.eps0; % move to dimension coinciding with El
%     El = permute(E{l}, [1 2 3 5 4]); % additional dimension for contraction with eps
%     el   = shiftdim( dat.gew.lay{l}.mat.e(udof,udof,udof), -3)/dat.gew.np.e0; % move to dimension for contraction
%     Fl = permute(F{l}, [1 2 3 6 4 5]);  % additional dimension for contraction with e
%     D{l} = sum(sum( el.*Fl ,5),6) + sum( epsl.*El ,5);
% end

end
