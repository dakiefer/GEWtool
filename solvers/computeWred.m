function w = computeWred(guw, ks)
    op = guw.op; % operators
    r = guw.geom.gdofBC(:); % degrees to remove
    k = setdiff([guw.geom.gdofOfLay{:}], r); % degrees to keep (unknowns)
    Mkk = op.M(k, k); 
    Mkr = op.M(k, r);
    L0kk = op.L0(k, k); L1kk = op.L1(k, k); L2kk = op.L2(k, k);
    L0kr = op.L0(k, r); L1kr = op.L1(k, r); L2kr = op.L2(k, r);
%     lay1 = guw.lay(1); lay2 = guw.lay(end);
%     opBl1 = lay1.tractionOp(0); opBl2 = lay2.tractionOp(0);
    B1k = op.L1(r,k); B0k = op.L0(r,k);
    B1r = op.L1(r,r); B0r = op.L0(r,r);

    tic 
    kh = ks*guw.np.h0;
    whn = nan(length(k), length(kh)); % change
    for ii = 1:length(kh)
        Bk = 1i*kh(ii)*B1k + B0k; Br = 1i*kh(ii)*B1r + B0r; 
        G = -Br\Bk;
        Lkk = (1i*kh(ii))^2*L2kk + (1i*kh(ii))*L1kk + L0kk;
        Lkr = (1i*kh(ii))^2*L2kr + (1i*kh(ii))*L1kr + L0kr;
        Lred = Lkk + Lkr*G;
        Mred = Mkk + Mkr*G;
        [wh2] = eig(Lred, -Mred); 
        whn(:,ii) = sqrt(wh2);
    end
%     whn(whn == 0) = nan;
    chron = toc;
    fprintf('nF: %d, nK: %d, elapsed time: %g, time per point: %g ms\n', size(whn, 2), size(whn, 1), chron, chron/length(whn(:))*1e3);
    w = whn*guw.np.fh0/guw.np.h0;
end
