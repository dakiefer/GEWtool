function [eelec] = energyDensityElectric(dat)
% energyDensityElectric - Compute the electric energy density.
% 
% eelec = 1/4*Re(E.D*)
%
% Usage: 
% > eelec = energyDensityElectric(dat);
% 
% Arguments: 
% - dat:    GEWdat object that stores the computed modes. 
%

if ~isscalar(dat) % compute recursively for every waveguide problem in the vector "dat"
    eelec = arrayfun(@energyDensityElectric,dat,'UniformOutput',false); % apply to every object in the array "dat"
    return;
end

E = electricField(dat); 
D = electricFluxDensity(dat); 

eelec = cell(dat.gew.geom.nLay,1);
for l = 1:dat.gew.geom.nLay
    eelec{l} = 1/4*real(sum(E{l}.*conj(D{l}), 4));
end

end
