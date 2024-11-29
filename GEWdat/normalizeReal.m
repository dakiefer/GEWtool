function datN = normalizeReal(dat)
% normalizeReal - Normalize propagating modes with real wavenumbers. 

if ~isscalar(dat) % compute recursively for every waveguide problem in the vector "dat"
    datN = arrayfun(@normalizeReal,dat); % apply to every object in the arrays "dat"
    return;
end

if any(dat.k ~= dat.k(:,1), 'all') % orthogonality for constant frequency! 
    error('GEWTOOL:normalizeReal:nonconstk', 'Use this normalization for modes at constant wavenumber only.');
end

P = powerFluxMag(dat); % power flux of each mode

% normalize each of the layers in the same way:
for l = 1:length(dat.gew.lay) 
    dat.u{l} = dat.u{l}./sqrt(abs(P)); % should be real anyways. ignore sign.
end

datN = dat;

end