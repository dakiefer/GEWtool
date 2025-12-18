function datRed = extractModes(dat, indk, indw)
% EXTRACTMODES - Returns the data of the indicated modes (wavenumber-frequency). 
% 
% Usage: 
% > datRed = extractModes(dat, indk, indw);
% or
% > datRed = extractModes(dat, sel);
%
% Arguments: 
% - dat:   dispersion data structure as returned by the solver functions
% - indk:  vector that indexes the desired wavenumbers (first dimension in dat fields)
% - indw:  vector that indexes the desired frequencies (second dimension in dat fields)
% - sel:   logical array of same size as dat.k or vector with linear indices.
%          Only entire rows or entire columns are removed if possible, the remaining
%          undesired data points are set to NaN.
% 
% Return value: 
% - datRed:   dispersion data structure consistent to dat but reduced to the
%             desired modes.
% 
% 2024-2025 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, CNRS, France

if ~isscalar(dat) % extract recursively for every waveguide problem in the vector "dat"
    extract = @(datObj) extractModes(datObj, indk, indw);
    datRed = arrayfun(extract,dat,'UniformOutput',true); % apply to every object in the array "dat"
    return;
end

if nargin < 3 
    if islogical(indk) && all(size(indk) == size(dat.k))
        sel = indk; 
    else % assume indk is a list of linear indices to keep.
        % We will actually keep entire rows and columns NaN remaining entries:
        sel = false(size(dat.k)); 
        sel(indk) = true; 
    end
    indw = any(sel,1);
    indk = any(sel,2);
    dat.k(~sel) = nan; 
    dat.w(~sel) = nan; 
end

datRed = GEWdat(dat.gew, dat.k(indk,indw), dat.w(indk,indw), dat.Psi(indk,indw,:)); 

end
