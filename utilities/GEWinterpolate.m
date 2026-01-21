function [ui, zi] = GEWinterpolate(gew, u, zi, parity)
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
% When ''gew'' is a symmetrized geometry with z in [0, zmax], GEWinterpolate
% will retrieve field data for zi < 0 by symmetrical extension. You just need to
% provide the negative coordinates. When providing the number of interpolation
% points instead of coordinates, GEWinterpolate will automatically interpolate
% on the full symmetric domain [-gew.geom.zItf(end), gew.geom.zItf(end)].
% 
% Example: 
% >> dat = computeW(gew,k);         % compute solutions for vector of wavenumber 'k'
% >> indk = 5;                      % fifth wavenumber
% >> indw = 6;                      % sixth frequency
% >> datMode = extractModes(dat,indk,indw); % data of only mode (indk, indw)
% >> [ui, zi] = GEWinterpolate(gew, datMode.u, 200); % interpolate onto 200 equi-distant coordinates zi
%
% Relies on the barycentric Lagrange interpolation implemented by Greg von Winckel 
% in ''barylag''.
%
% See also: extractModes, barylag
%
% 2024-2026 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, CNRS, France

if ~iscell(u) && gew.geom.nLay ~= 1  % for multilayered waveguide we need a cell array for u
    error('GEWTOOL:GEWinterpolate:incorrectDataStructure', 'The waveguide has multiple layers but the data structure is not a cell array.');
elseif ~iscell(u) && gew.geom.nLay == 1
    uu{1} = u; u = uu; % wrap into cell for consistency
end
if nargin >= 4 && ~(parity == 1 || parity == -1)
    error('GEWTOOL:GEWinterpolate','The forth arguement "parity" should be 1 or -1.')
end

% convert from number of nodes Ni to nodes zi, if necessary:
if isscalar(zi) && isreal(zi) && mod(zi,1) == 0 % test if integer real scalar
    Ni = zi; % the number of nodes was provided rather than the interpolation grid.
    if gew.geom.symmetrized 
        zi = linspace(-gew.geom.zItf(end), gew.geom.zItf(end), Ni).';
    else
        zi = linspace(gew.geom.zItf(1), gew.geom.zItf(end), Ni).';
    end
end

% ensure correct structure of data:
for l = 1:gew.geom.nLay
    u{l} = squeeze(u{l});  % extract layer data and squeeze
    N = gew.geom.N(l);     % expected number of nodes 
    Ndata = size(u{l},1);  % actual number of nodes
    if N ~= Ndata 
        error('GEWTOOL:GEWinterpolate',...
        'Size inconsistency. You must reduce the data to one mode before interpolating. You can use extractModes() for this purpose. See "help extractModes".');
    end
end

% initialize:
zi = zi(:); % barylag needs a column vector
s = size(u{1}); s(1) = length(zi); % size of data structure after interpolation
ui = zeros(s); % allocate
nComp = prod(s(2:end)); % number of components of the field (3 for displ, 9 for stress...)

% % guess parity of functions: 
if nargin < 4 && gew.geom.symmetrized 
    if gew.geom.gdofDBC(1) == 1 % anti-symmetric waves
        parity = +1;  % this depends on the phase convention of eig vecs in Matlab
    else % symmetric waves
        parity = -1;  % this depends on the phase convention of eig vecs in Matlab
    end
    if nComp >= 4 % strain and stress have the opposit parity to the displacements
        parity = -1*parity; 
    end
end

% interpolate+extrapolate on every layer and every field component:
for l = 1:gew.geom.nLay 
    for n = 1:nComp % loop over all component of the data 
        zl = gew.geom.z{l};  % nodal points
        datal = [zl, u{l}(:,n)];               % initial data
        indl = zi >= zl(1) & zi <= zl(end);    % indices of interpolated data for this layer
        ui(indl,n) = barylag(datal, zi(indl)); % interpolate onto corresponding zi
        if gew.geom.symmetrized && any(zi < 0) % extend symmetrically to negative coordinates
            indb = zi >= -zl(end) & zi <= -zl(1); % indices of interpolation points on symmetric range [-zl(1) -zl(end)]
            if all(-flip(zi(indb)) == zi(indl)) % interpolation points on positive and negative range are the same: we can use the last interpolation
                indRev = flip(find(indl));
                ui(indb,n) =  parity*conj(ui(indRev,n));
            else % we need to interpolate the negative section separately
                uTmp = barylag(datal, -flip(zi(indb))); % interpolate onto corresponding zi
                ui(indb,n) =  parity*conj(flip(uTmp));
            end
        end
    end
end

end
