function dat = computeW(guw, k, nModes)
    if ~isvector(k), error('Wavenumbers should be a [1xN] array.'); end
    if nargin < 3, nModes = size(guw.op.M,1); end
    k = k(:).'; % row vector
    kh = k*guw.np.h0;
    M = guw.op.M; L0 = guw.op.L0; L1 = guw.op.L1; L2 = guw.op.L2;
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
    dat.w = whn*guw.np.fh0/guw.np.h0;
    dat.k = kh.*ones(size(whn))/guw.np.h0;
    dat.u = cell(guw.geom.nLay, 1);
    for i = 1:guw.geom.nLay
        ulay = u(:,:,guw.geom.gdofOfLay{i});
        dat.u{i} = reshape(ulay, [size(whn), guw.geom.N(i), guw.geom.Nudof(i)]);
    end
end
