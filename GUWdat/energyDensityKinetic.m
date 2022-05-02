function [ekin] = energyDensityKinetic(wguide, dat)
%ENERGYDENSITYKINETIK Compute the kinetik energy density Ekin.

rho = wguide.lay.mat.rho;
v = velocity(dat); 
ekin = 1/4*rho*real(sum(conj(v).*v, 4));

end
