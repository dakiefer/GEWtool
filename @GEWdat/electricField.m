function [E] = electricField(dat)
% electricField - returns the electric field intensity E. 
% 
% E = -grad(phi) = - (ex ik phi + ez ∂z phi)
% 
% Only relevant for piezoelectric waveguides. A zero-cell array of appropriate size 
% is returned otherwise. 
% 
% Usage: 
% > E = electricField(dat)
% 
% Arguments: 
% - dat:    GEWdat object that stores the computed modes. 
%

if ~isscalar(dat) % compute recursively for every waveguide problem in the vector "gew"
    E = arrayfun(@electricField,dat,'UniformOutput',false); % apply to every object in the array "dat"
    return;
end

G = potentialGrad(dat.gew,dat); 

E = cell(dat.gew.geom.nLay, 1); % allocate for each layer
for l = 1:dat.gew.geom.nLay
    E{l} = -G{l}; 
end

end
