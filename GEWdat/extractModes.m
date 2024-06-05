function datRed = extractModes(dat, indk, indw)
% EXTRACTMODES - Returns the data of the indicated modes (wavenumber-frequency). 
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

nLay = length(dat.u); % number of layers
uRed = cell(1,nLay); % allocate
for l = 1:nLay % for every layer
    uRed{l} = dat.u{l}(indk,indw,:,:);
end
datRed.w = dat.w(indk,indw); 
datRed.k = dat.k(indk,indw);
datRed.u = uRed; 

end
