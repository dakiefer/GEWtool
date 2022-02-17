function [w, u] = computeWandU(guide, k, nModes)
    if nargin < 3, nModes = size(guide.op.M,1); end
    kh = k*guide.np.h0;
    M = guide.op.M; L0 = guide.op.L0; L1 = guide.op.L1; L2 = guide.op.L2;
    tic 
    whn = nan(size(M,1), length(kh));
    u = nan(nModes, length(kh), size(M,1));
    for ii = 1:length(kh)
        [ui, wh2] = polyeig((1i*kh(ii))^2*L2 + (1i*kh(ii))*L1 + L0, M); % does not work properly with eig()
        [whi, ind] = sort(real(sqrt(wh2)));
        ui = ui(:,ind);
        whn(:,ii) = whi(1:nModes);
        u(:,ii,:) = (ui(:,1:nModes)).';
    end
    whn(whn == 0) = nan;
    chron = toc;
    fprintf('nF: %d, nK: %d, elapsed time: %g, time per point: %g ms\n', size(whn, 2), size(whn, 1), chron, chron/length(whn(:))*1e3);
    w = whn*guide.np.fh0/guide.np.h0;
end
