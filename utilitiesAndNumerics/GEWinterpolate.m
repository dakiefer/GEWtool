function [ui, yi] = GEWinterpolate(gew, u, yi)
% GEWinterpolate - Interpolate the data over all layers of the waveguide.
% 
% Arguments: 
% - gew:     (Waveguide object) Description of the waveguide geometry and
%            discretization
% - u:       (cell array) Data of _one_ selected mode to be interpolated,
%            e.g., ''dat.u(1,1,:,:)'' returned by the solvers.
% - yi:      - scalar: numer of coordinates that will be generated on an
%              equidistant grid on [gew.geom.yItf(1), gew.geom.yItf(end)].
%            - vector: the coordinates on which to interpolate (should be 
%              in the domain of the waveguide cross-section.
% 
% When ''gew'' is a symmetrized geometry, i.e., y in [0, ymax], GEWinterpolate
% will retrieve field data for yi < 0 by symmetrical extension. You just need to
% provide the negative coordinates. When providing the number of interpolation
% points instead of coordinates, GEWinterpolate will automatically interpolate
% on the symmetric domain [-gew.geom.yItf(end), gew.geom.yItf(end)].
%
% Relies on the barycentric Lagrange interpolation implemented by Greg von
% Winckel in ''barylag''.
%
% See also: barylag
%
% 2024 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France

if isscalar(yi) && isreal(yi) && mod(yi,1) == 0 % test if integer real scalar
    Ni = yi; % the number of nodes was provided rather than the interpolation grid.
    if gew.geom.symmetrized 
        yi = linspace(-gew.geom.yItf(end), gew.geom.yItf(end), Ni).';
    else
        yi = linspace(gew.geom.yItf(1), gew.geom.yItf(end), Ni).';
    end
end

for l = 1:gew.geom.nLay
    u{l} = squeeze(u{l});  % remove singleton dimension remaining from selecting a mode
    N = gew.geom.N(l);     % expected number of nodes 
    Ndata = size(u{l},1);  % actual number of nodes
    if N ~= Ndata 
        error('GEWTOOL:GEWinterpolate',...
        'Size inconsistency. You must reduce the data to one mode before interpolating. You can use extractModes() for this purpose. See "help extractModes".');
    end
end

yi = yi(:); % barylag needs a column vector
s = size(u{1}); s(1) = length(yi); % size of data structure after interpolation
ui = zeros(s); % allocate
for l = 1:gew.geom.nLay 
    for n = 1:prod(s(2:end)) % loop over all component of the data 
        yl = gew.geom.y{l};  % nodal points
        datal = [yl, u{l}(:,n)];               % initial data
        indl = yi >= yl(1) & yi <= yl(end);    % indices of interpolated data for this layer
        ui(indl,n) = barylag(datal, yi(indl)); % interpolate onto corresponding yi
        if gew.geom.symmetrized && any(yi < 0) % extend symmetrically to negative coordinates
            indb = yi >= -yl(end) & yi <= -yl(1);
            if all(yi(indb) == yi(indl)) % we can use the last interpolation
                indRev = flip(find(indl));
                if gew.geom.gdofDBC == gew.geom.N+1 % symmetric waves
                    ui(indb,n) = -conj(ui(indRev,n));
                else % anti-symmetric waves
                    ui(indb,n) = conj(ui(indRev,n));
                end
            else % we need to interpolate the negative section separately
                uTmp = barylag(datal, -flip(yi(indb))); % interpolate onto corresponding yi
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
