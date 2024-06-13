function dat = normalizeComplex(gew, dat)
% normalizeToPowerFlux Normalize the displacement eigenvectors to unit power flux
% of the waves. Complex waves are normalized by the cross power flux with their
% complex conjugate pair. 
%
% B. A. Auld, Acoustic Fields and Waves in Solids 2, 2nd ed., vol. 2, 2 vols. 
% Malabar: Krieger Publishing Company, 1990.

warning('GEWTOOL:untested', 'Untested function.')

if any(dat.w ~= dat.w(:,1), 'all') % orthogonality for constant frequency! 
    error('GEWTOOL:normalizeComplex:nonconstfreq', 'Cross power flux must be computed for wavenumbers at constant frequency.');
end

v = velocity(dat);
T = stress(gew, dat);

for m = 1:size(dat.k, 2)
    % search "sibling" mode indkn corresponding to indkm:
    indm = sub2ind(size(dat.k), 1:size(dat.k,2), m*ones(1,size(dat.k,2)));
    km = dat.k(indm);
    [~, indn] = min(abs(dat.k - conj(km.')), [], 1); % TODO: test if this works for multiple frequencies
    indn = sub2ind(size(dat.k), 1:size(dat.k,2), indn.');
    % compute the power flux Pmn between these "sibling" modes (do not used
    % crossPowerFlux(), as this would calculate Pmn for all m and n):
    Imn = cell(1, gew.geom.nLay);
    for l = 1:length(gew.lay)
%         indkn = [indkn, (numel(dat.k)+1):numel(v{1})];
        vLay = v{l}; vLay = reshape(vLay, numel(dat.k), gew.geom.N(l), gew.geom.Nudof(l));
        TLay = T{l}; TLay = reshape(TLay, numel(dat.k), gew.geom.N(l), gew.geom.Nudof(l), gew.geom.Nudof(l));
        vm = vLay(:,indm,:);
        vn = vLay(:,indn,:);
        txm = permute(TLay(:,indm,1,:), [1 2 4 3]); % traction ex.T = T'.ex of size [nF x nK x yi x tx]
        txn = permute(TLay(:,indn,1,:), [1 2 4 3]); % traction ex.T = T'.ex of size [nF x nK x yi x tx]
        Imn{l} = sum(-conj(vn).*txm - vm.*conj(txn), 3); % power flux density
    end
    Pmn = 1/4*GEWintegrate(gew, Imn, [], 2); % size: [nF]
    % normalize mode ik:
    for l = 1:length(gew.lay) 
        um = dat.u{l}(m,:,:,:);
        dat.u{l}(m,:,:,:) = um./sqrt(abs(Pmn)); % should be real anyways. ignore sign.
    end
end

end
