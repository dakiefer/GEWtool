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
    ekin = arrayfun(@energyDensityKinetic,dat,'UniformOutput',false); % apply to every object in the array "dat"
    return;
end

v = velocity(dat); 

ekin = cell(dat.gew.geom.nLay,1);
for l = 1:dat.gew.geom.nLay
    rho = dat.gew.lay{l}.mat.rho/dat.gew.np.rho0;
    ekin{l} = 1/4*rho*real(sum(conj(v{l}).*v{l}, 4));
end

end
