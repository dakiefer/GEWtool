function dat = computeK(wguide, w, nModes)
    if nargin < 3, nModes = 2*size(wguide.op.M,1); end
    wh = w*wguide.np.h0;
    M = wguide.op.M; L0 = wguide.op.L0; L1 = wguide.op.L1; L2 = wguide.op.L2;
    tic 
    kh = nan(length(wh), nModes);
    u = nan(length(wh), nModes, size(M,1));
    for ii = 1:length(wh)
        whn = wh(ii)/wguide.np.fh0; % current frequency-thickness (normalized)
        [ui, ikhi] = polyeig(L0 + whn^2*M, L1, L2);
        spurious = abs(ikhi)<=1e-8; ikhi(spurious) = nan; ui(:,spurious) = nan;
        [khi, ind] = sort(-1i*ikhi);
        ui = ui(:,ind); % sort
        u(ii,:,:) = ui(:,1:nModes).'; % save
        kh(ii, :) = khi(1:nModes);
    end
    chron = toc; 
    fprintf('nF: %d, nK: %d, elapsed time: %g, time per point: %g. ms\n', size(kh, 1), size(kh, 2), chron, chron/length(kh(:))*1e3);
    dat.k = kh/wguide.np.h0;
    dat.u = reshape(u, [size(kh), wguide.geom.N, wguide.geom.Nudof]); 
    dat.w = w.*ones(size(kh));
end
