function [ui, zi] = GEWinterpolate(gew, u, zi)
% GEWinterpolate - Interpolate the data over all layers of the waveguide.
% 
% Arguments: 
% - gew:      (Waveguide object) Description of the waveguide.
% - u:        (cell array of numeric arrays | numeric array) Data of _one_
%             selected mode to be interpolated, e.g., ''dat.u{1}(1,1,:,:)'' returned 
%             by computeW() or computeK(). Use extractModes() to extract the data of 
%             one mode of a mulit-layer waveguide.
% - Ni | zi:  - Ni (scalar): numer of coordinates that will be generated on an
%               equidistant grid on [gew.geom.zItf(1), gew.geom.zItf(end)].
%             - zi (vector): the coordinates on which to interpolate (should be 
%               in the domain of the provided data (i.e. the waveguide cross-section).
% 
% Return values: 
% - ui:    (array) data interpolated to coordinates zi
% - zi:    (vector) coordinates zi
% 
% When ''gew'' is a symmetrized geometry, i.e., z in [0, zmax], GEWinterpolate
% will retrieve field data for zi < 0 by symmetrical extension. You just need to
% provide the negative coordinates. When providing the number of interpolation
% points instead of coordinates, GEWinterpolate will automatically interpolate
% on the symmetric domain [-gew.geom.zItf(end), gew.geom.zItf(end)].
% 
% Example: 
% >> dat = computeW(gew,k);         % compute solutions for vector of wavenumber 'k'
% >> indk = 5;                      % fifth wavenumber
% >> indw = 6;                      % sixth frequency
% >> datMode = extractModes(dat,indk,indw); % data of only mode (indk, indw)
% >> [ui, zi] = GEWinterpolate(gew, datMode.u, 200); % interpolate onto 200 equi-distant coordinates zi
%
% Relies on the barycentric Lagrange interpolation implemented by Greg von
% Winckel in ''barylag''.
%
% See also: extractModes, barylag
%
% 2024 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France

if ~iscell(u) && gew.geom.nLay ~= 1  % for multilayered waveguide we need a cell array for u
    error('GEWTOOL:GEWinterpolate:incorrectDataStructure', 'The waveguide has multiple layers but the data structure is not a cell array.');
elseif ~iscell(u) && gew.geom.nLay == 1
    uu{1} = u; u = uu; % wrap into cell for consistency
end
if isscalar(zi) && isreal(zi) && mod(zi,1) == 0 % test if integer real scalar
    Ni = zi; % the number of nodes was provided rather than the interpolation grid.
    if gew.geom.symmetrized 
        zi = linspace(-gew.geom.zItf(end), gew.geom.zItf(end), Ni).';
    else
        zi = linspace(gew.geom.zItf(1), gew.geom.zItf(end), Ni).';
    end
end

for l = 1:gew.geom.nLay
    u{l} = squeeze(u{l});  % extract layer data and squeeze
    N = gew.geom.N(l);     % expected number of nodes 
    Ndata = size(u{l},1);  % actual number of nodes
    if N ~= Ndata 
        error('GEWTOOL:GEWinterpolate',...
        'Size inconsistency. You must reduce the data to one mode before interpolating. You can use extractModes() for this purpose. See "help extractModes".');
    end
end

zi = zi(:); % barylag needs a column vector
s = size(u{1}); s(1) = length(zi); % size of data structure after interpolation
ui = zeros(s); % allocate
for l = 1:gew.geom.nLay 
    for n = 1:prod(s(2:end)) % loop over all component of the data 
        zl = gew.geom.z{l};  % nodal points
        datal = [zl, u{l}(:,n)];               % initial data
        indl = zi >= zl(1) & zi <= zl(end);    % indices of interpolated data for this layer
        ui(indl,n) = barylag(datal, zi(indl)); % interpolate onto corresponding zi
        if gew.geom.symmetrized && any(zi < 0) % extend symmetrically to negative coordinates
            indb = zi >= -zl(end) & zi <= -zl(1);
            if all(zi(indb) == zi(indl)) % we can use the last interpolation
                indRev = flip(find(indl));
                if gew.geom.gdofDBC == gew.geom.N+1 % symmetric waves
                    ui(indb,n) = -conj(ui(indRev,n));
                else % anti-symmetric waves
                    ui(indb,n) = conj(ui(indRev,n));
                end
            else % we need to interpolate the negative section separately
                uTmp = barylag(datal, -flip(zi(indb))); % interpolate onto corresponding zi
                if gew.geom.gdofDBC == gew.geom.N+1 % symmetric waves
                    ui(indb,n) = -conj(flip(uTmp));
                else % anti-symmetric waves
                    ui(indb,n) = conj(flip(uTmp));
                end
            end
        end
    end
end

end
