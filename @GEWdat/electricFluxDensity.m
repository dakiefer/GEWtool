function [D] = electricFluxDensity(dat)
% electricFluxDensity - returns the electric flux density D. 
% 
% D = eps*E 
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

E = electricField(dat); 
udof = dat.gew.udof; % E has same components as displacements u

D = cell(dat.gew.geom.nLay, 1); % allocate for each layer
for l = 1:dat.gew.geom.nLay
    epsl = dat.gew.lay{l}.mat.epsilon(udof,udof);
    El = permute(E{l}, [1 2 3 5 4]);         % additional dimension for contraction with eps
    eepsl = shiftdim(epsl, -3);   % move to dimension coinciding with El
    D{l} = sum(eepsl.*El,5); 
end

end
