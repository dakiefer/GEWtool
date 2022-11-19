function dat = computeK(gews, w, nModes)
    % computeK - Obtain complex wavenumbers k for specified frequencies w.
    % Solves the polynomial eigenvalue problem [(ik)^2*L2 + ik*L1 + L0(w)]*u = 0.
    %
    % Arguments:
    % - gews:    Waveguide object(s), either a scalar or vector.
    %            Describes the eigenproblem, i.e., the matrices Li.
    %            If gews is a vector, computeW solves one problem after another 
    %            and returns a vector of results "dat" of same length.
    % - w:       Angular frequencies to specify in rad/s. Vector valued.
    % - nModes:  (optional) Number of modes to save (discards the highest wavenumbers).
    %
    % Return value:
    % - dat:     A data structure containing 
    %            - w: the angular frequencies in rad/s, expanded to [nF x nK]
    %            - k: the wavenumbers in rad/m [nF x nK]
    %            - u: the displacement eigenvectors as a 
    %                 cell array describing the layers, elements are [nF x nK x N x Nudof]
    % 
    % See also computeW, Waveguide.
    % 
    % 2022 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France
    
    if ~isvector(w), error('Angular frequencies should be a [Nx1] array.'); end
    w = w(:); % column vector
    for i = 1:length(gews) % solve for a list of waveguide objects
        gew = gews(i);
        wh = w*gew.np.h0;
        M = gew.op.M; L0 = gew.op.L0; L1 = gew.op.L1; L2 = gew.op.L2;
        if nargin < 3
            nModes = size(gew.op.M,1);  % for now we dont distinguish between linearized and quadratic EVP
        end
        if nModes > size(gew.op.M,1)
            warning('GEWTOOL:computeK:tooManyModes', 'More modes requested than available. Resetting nModes to the matrix size.')
            nModes = size(gew.op.M,1);
        end
        kh = nan(length(wh), nModes);
        u = nan(length(wh), nModes, gew.geom.Ndof);
        gdoffree = setdiff([gew.geom.gdofOfLay{:}], gew.geom.gdofDBC(:).');
        for ii = 1:length(wh)
            whn = wh(ii)/gew.np.fh0; % current frequency-thickness (normalized)
            if isequal(L2, zeros(size(L2))) % is linearized as [(ik) L1 + L0 + w^2 M]*u = 0
                [ui, ikhi] = polyeig(L0 + whn^2*M, L1); 
            elseif isempty(L1) % is linearized as [(ik)^2 L2 + L0 + w^2 M]*u = 0
                [ui, ikhi2] = polyeig(L0 + whn^2*M, L2); % calculate (ikh)^2
                ikhi = sqrt(ikhi2);
                N = gew.geom.N; dofy = N+1:2*N;
                ui(dofy,:) = ikhi.'.*ui(dofy,:); % eig.vec. was [ux, 1i*k*uy]
            else % quadratic EVP: [(ik)^2 L2 + (ik) L1 + L0 + w^2 M]*u = 0
                [ui, ikhi] = polyeig(L0 + whn^2*M, L1, L2);
            end
            spurious = isinf(ikhi) | isnan(ikhi); ikhi(spurious) = []; ui(:,spurious) = [];
            [khi, ind] = sort(-1i*ikhi);
            ui = ui(:,ind); % sort
    %         [~, ind] = sort(real(khi));
    %         khi = khi(ind); ui = ui(:,ind); % sort
            u(ii,:,gdoffree) = ui(:,1:nModes).'; % save
            kh(ii, :) = khi(1:nModes);
        end
        dat(i).k = kh/gew.np.h0;
        dat(i).w = w.*ones(size(kh));
        dat(i).u = cell(gew.geom.nLay, 1); % initialize
        for l = 1:gew.geom.nLay
            ulay = u(:,:,gew.geom.gdofOfLay{l});
            dat(i).u{l} = reshape(ulay, [size(kh), gew.geom.N(l), gew.geom.Nudof(l)]);
        end
    end
end
