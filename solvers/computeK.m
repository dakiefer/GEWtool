function dat = computeK(guw, w, nModes)
    if ~isvector(w), error('Angular frequencies should be a [Nx1] array.'); end
    w = w(:); % column vector
    wh = w*guw.np.h0;
    M = guw.op.M; L0 = guw.op.L0; L1 = guw.op.L1; L2 = guw.op.L2;
    if nargin < 3
        if isequal(L2, zeros(size(L2)))
            nModes = size(guw.op.M,1);
        else
            nModes = 2*size(guw.op.M,1);
        end
    end
    tic 
    kh = nan(length(wh), nModes);
    u = nan(length(wh), nModes, size(M,1));
    for ii = 1:length(wh)
        whn = wh(ii)/guw.np.fh0; % current frequency-thickness (normalized)
        if isequal(L2, zeros(size(L2)))
            [ui, ikhi] = polyeig(L0 + whn^2*M, L1);
        else
            [ui, ikhi] = polyeig(L0 + whn^2*M, L1, L2);
        end
        spurious = abs(ikhi)<=1e-8; ikhi(spurious) = nan; ui(:,spurious) = nan;
        [khi, ind] = sort(-1i*ikhi);
        ui = ui(:,ind); % sort
        u(ii,:,:) = ui(:,1:nModes).'; % save
        kh(ii, :) = khi(1:nModes);
    end
    chron = toc; 
    fprintf('nF: %d, nK: %d, elapsed time: %g, time per point: %g. ms\n', size(kh, 1), size(kh, 2), chron, chron/length(kh(:))*1e3);
    dat.k = kh/guw.np.h0;
    dat.w = w.*ones(size(kh));
    dat.u = cell(guw.geom.nLay, 1);
    for i = 1:guw.geom.nLay
        ulay = u(:,:,guw.geom.gdofOfLay{i});
        dat.u{i} = reshape(ulay, [size(kh), guw.geom.N(i), guw.geom.Nudof(i)]);
    end
end
