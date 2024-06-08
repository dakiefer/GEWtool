function datN = normalizeReal(gew,dat)
% normalizeReal - Normalize propagating modes with real wavenumbers. 

if any(dat.k ~= dat.k(:,1), 'all') % orthogonality for constant frequency! 
    error('GEWTOOL:normalizeReal:nonconstk', 'Use this normalization for constant constant wavenumber computations only.');
end

P = powerFluxMag(gew, dat); % power flux of each mode

% normalize each of the layers in the same way:
for l = 1:length(gew.lay) 
    dat.u{l} = dat.u{l}./sqrt(abs(P)); % should be real anyways. ignore sign.
end

datN = dat;

end