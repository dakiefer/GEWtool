function datRed = extractModes(dat, indk, indw)
% EXTRACTMODES - Returns the data of the indicated modes (wavenumber-frequency). 
% 
% Usage: 
% > datRed = extractModes(dat, indk, indw);
%
% Arguments: 
% - dat:   dispersion data structure as returned by the solver functions
% - indk:  vector that indexes the desired wavenumbers (first dimension in dat fields)
% - indw:  vector that indexes the desired frequencies (second dimension in dat fields)
% 
% Return value: 
% - datRed:   dispersion data structure consistent to dat but reduced to the
%             desired modes.
% 
% 2024 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France

if ~isscalar(dat) % compute recursively for every waveguide problem in the vector "dat"
    compute = @(datObj) extractModes(datObj, indk, indw);
    datRed = arrayfun(compute,dat,'UniformOutput',true); % apply to every object in the array "dat"
    return;
end

datRed = GEWdat(dat.gew, dat.k(indk,indw), dat.w(indk,indw), dat.Psi(indk,indw,:)); 

end
