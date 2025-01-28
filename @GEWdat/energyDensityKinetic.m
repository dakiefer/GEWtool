function [ekin] = energyDensityKinetic(dat)
%ENERGYDENSITYKINETIK Compute the kinetik energy density Ekin.
%
% Usage: 
% > ekin = energyDensityKinetic(dat);
% 
% Arguments: 
% - dat:    GEWdat object that stores the computed modes. 
%

if ~isscalar(dat) % compute recursively for every waveguide problem in the vector "dat"
    ekin = arrayfun(@energyDensityKinetic,dat,'UniformOutput',false); % apply to every object in the arrays "dat"
    return;
end

v = velocity(dat); 

ekin = cell(dat.gew.geom.nLay,1);
for i = 1:dat.gew.geom.nLay
    rho = dat.gew.lay{i}.mat.rho;
    ekin{i} = 1/4*rho*real(sum(conj(v{i}).*v{i}, 4));
end

end
