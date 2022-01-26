function whn = solveWithKh(kh, M, L0, L1, L2)
    tic 
    whn = nan(size(M, 2), length(kh));
    for ii = 1:length(kh)
        [wh2] = polyeig((1i*kh(ii))^2*L2 + (1i*kh(ii))*L1 + L0, M); % does not work properly with eig()
        whn(:,ii) = real(sqrt(wh2));
    end
    whn(whn == 0) = nan;
    chron = toc;
    fprintf('nF: %d, nK: %d, elapsed time: %g, time per point: %g. ms\n', size(whn, 2), size(whn, 1), chron, chron/length(whn(:))*1e3);
end
