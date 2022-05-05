function [ekin] = energyDensityKinetic(wguide, dat)
%ENERGYDENSITYKINETIK Compute the kinetik energy density Ekin.

v = velocity(dat); 

ekin = cell(wguide.geom.nLay,1);
for i = 1:wguide.geom.nLay
    rho = wguide.lay(i).mat.rho;
    ekin{i} = 1/4*rho*real(sum(conj(v{i}).*v{i}, 4));
end

end
