function dat = normalizeToPowerFlux(guw, dat)
% normalizeToPowerFlux Normalize the displacement eigenvectors to unit power flux
% of the waves. Complex waves are normalized by the cross power flux with their
% complex conjugate pair. 
%
% B. A. Auld, Acoustic Fields and Waves in Solids 2, 2nd ed., vol. 2, 2 vols. 
% Malabar: Krieger Publishing Company, 1990.

if any(dat.w ~= dat.w(:,1), 'all') % orthogonality for constant frequency! 
    error('GEWTOOL:crossPowerFlux:nonconstfreq', 'Cross power flux must be computed for wavenumbers at constant frequency.');
end

v = velocity(dat);
T = stress(guw, dat);

for ik = 1:size(dat.k, 2)
    indkm = sub2ind(size(dat.k), 1:size(dat.k,1), ik*ones(1,size(dat.k,1)));
    km = dat.k(indkm);
    [~, indkn] = min(abs(dat.k - conj(km.')), [], 2); % TODO: test if this works for multiple frequencies
    indkn = sub2ind(size(dat.k), 1:size(dat.k,1), indkn.');
    Imn = cell(1, guw.geom.nLay);
    for l = 1:length(guw.lay)
%         indkn = [indkn, (numel(dat.k)+1):numel(v{1})];
        vLay = v{l}; vLay = reshape(vLay, numel(dat.k), guw.geom.N(l), guw.geom.Nudof(l));
        TLay = T{l}; TLay = reshape(TLay, numel(dat.k), guw.geom.N(l), guw.geom.Nudof(l), guw.geom.Nudof(l));
        vm = vLay(indkm,:,:);
        vn = vLay(indkn,:,:);
        txm = permute(TLay(indkm,:,1,:), [1 2 4 3]); % traction ex.T = T'.ex of size [nF x nK x yi x tx]
        txn = permute(TLay(indkn,:,1,:), [1 2 4 3]); % traction ex.T = T'.ex of size [nF x nK x yi x tx]
        Imn{l} = sum(-conj(vn).*txm - vm.*conj(txn), 3); % power flux density
    end
    Pmn = 1/4*GEWintegrate(guw, Imn, 2); % size: [nF]
    % normalize:
    for l = 1:length(guw.lay) 
        um = dat.u{l}(:,ik,:,:);
        dat.u{l}(:,ik,:,:) = um./sqrt(abs(real(Pmn))); % should be real anyways. ignore sign.
    end
end

end
