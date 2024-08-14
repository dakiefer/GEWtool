function [ekin] = energyDensityKinetic(gew, dat)
%ENERGYDENSITYKINETIK Compute the kinetik energy density Ekin.

if ~isscalar(gew) % compute recursively for every waveguide problem in the vector "gew"
    compute = @(gewObj,datObj) energyDensityKinetic(gewObj, datObj); % function to apply
    ekin = arrayfun(compute,gew,dat,'UniformOutput',false); % apply to every object in the arrays "gew" and "dat"
    return;
end

v = velocity(dat); 

ekin = cell(gew.geom.nLay,1);
for i = 1:gew.geom.nLay
    rho = gew.lay{i}.mat.rho;
    ekin{i} = 1/4*rho*real(sum(conj(v{i}).*v{i}, 4));
end

end
