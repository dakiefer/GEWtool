function dat = computeW(wguide, k, nModes)
    if nargin < 3, nModes = size(wguide.op.M,1); end
    kh = k*wguide.np.h0;
    M = wguide.op.M; L0 = wguide.op.L0; L1 = wguide.op.L1; L2 = wguide.op.L2;
    tic 
    whn = nan(size(M, 2), length(kh));
    u = nan(nModes, length(kh), size(M,1));
    for n = 1:length(kh)
        [un, whn2] = polyeig((1i*kh(n))^2*L2 + (1i*kh(n))*L1 + L0, M); % does not work properly with eig()
        whnn = real(sqrt(whn2));
        spurious = whnn==0; whnn(spurious) = nan; un(:,spurious) = nan;
        [whnn, ind] = sort(whnn);
        un = un(:,ind);
        u(:,n,:) = un(:,1:nModes).'; % save
        whn(:, n) = whnn(1:nModes);
    end
    chron = toc;
    fprintf('nF: %d, nK: %d, elapsed time: %g, time per point: %g ms\n', size(whn, 2), size(whn, 1), chron, chron/length(whn(:))*1e3);
    dat.w = whn*wguide.np.fh0/wguide.np.h0;
    dat.k = kh.*ones(size(whn))/wguide.np.h0;
%     dat.u = reshape(u, [size(whn), wguide.geom.N, wguide.geom.Nudof]);
end
