function [ekin] = energyDensityKinetic(gew, dat)
%ENERGYDENSITYKINETIK Compute the kinetik energy density Ekin.

v = velocity(dat); 

ekin = cell(gew.geom.nLay,1);
for i = 1:gew.geom.nLay
    rho = gew.lay(i).mat.rho;
    ekin{i} = 1/4*rho*real(sum(conj(v{i}).*v{i}, 4));
end

end
