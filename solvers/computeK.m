function dat = computeK(guw, w, nModes)
    if ~isvector(w), error('Angular frequencies should be a [Nx1] array.'); end
    w = w(:); % column vector
    wh = w*guw.np.h0;
    M = guw.op.M; L0 = guw.op.L0; L1 = guw.op.L1; L2 = guw.op.L2;
    if nargin < 3
        nModes = size(guw.op.M,1);  % for now we dont distinguish between linearized and quadratic EVP
    end
    kh = nan(length(wh), nModes);
    u = nan(length(wh), nModes, size(M,1));
    for ii = 1:length(wh)
        whn = wh(ii)/guw.np.fh0; % current frequency-thickness (normalized)
        if isequal(L2, zeros(size(L2))) % is linearized as [(ik) L1 + L0 + w^2 M]*u = 0
            [ui, ikhi] = polyeig(L0 + whn^2*M, L1); 
        elseif isempty(L1) % is linearized as [(ik)^2 L2 + L0 + w^2 M]*u = 0
            [ui, ikhi2] = polyeig(L0 + whn^2*M, L2); % calculate (ikh)^2
            ikhi = sqrt(ikhi2);
            N = guw.geom.N; dofy = N+1:2*N;
            ui(dofy,:) = ikhi.'.*ui(dofy,:); % eig.vec. was [ux, 1i*k*uy]
        else % quadratic EVP: [(ik)^2 L2 + (ik) L1 + L0 + w^2 M]*u = 0
            [ui, ikhi] = polyeig(L0 + whn^2*M, L1, L2);
        end
        spurious = isinf(ikhi) | isnan(ikhi); ikhi(spurious) = []; ui(:,spurious) = [];
        [khi, ind] = sort(-1i*ikhi);
        ui = ui(:,ind); % sort
%         [~, ind] = sort(real(khi));
%         khi = khi(ind); ui = ui(:,ind); % sort
        u(ii,:,:) = ui(:,1:nModes).'; % save
        kh(ii, :) = khi(1:nModes);
    end
    dat.k = kh/guw.np.h0;
    dat.w = w.*ones(size(kh));
    dat.u = cell(guw.geom.nLay, 1); % initialize
    for i = 1:guw.geom.nLay
        ulay = u(:,:,guw.geom.gdofOfLay{i});
        dat.u{i} = reshape(ulay, [size(kh), guw.geom.N(i), guw.geom.Nudof(i)]);
    end
end
