function [ce] = energyVel(wguide, dat)
%ENERGYVEL Compute the energy velocity ce.

P = powerFlux(wguide, dat);
H = energyTotal(wguide, dat);
ce = P./H;

end
