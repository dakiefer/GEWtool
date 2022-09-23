function dat = computeW(guws, k, nModes)
    % computeW - Obtain frequency w for specified wavenumbers k.
    % Solves the eigenvalue problem described by L(k)*u = -w^2*M*u.
    %
    % Arguments:
    % - guws:    Waveguide object(s), either a scalar or vector.
    %            Describes the eigenproblem, i.e., the matrices L and M.
    %            If guws is a vector, computeW solves one problem after another 
    %            and returns a vector of results "dat" of same length.
    % - k:       Wavenumbers to specify in rad/m. Vector valued.
    % - nModes:  (optional) Number of modes to save (discards the highest frequencies).
    %
    % Return value:
    % - dat:     A data structure containing 
    %            - w: the angular frequencies in rad/s [nF x nK]
    %            - k: the wavenumbers in rad/m, expanded to [nF x nK]
    %            - u: the displacement eigenvectors as a 
    %                 cell array describing the layers, elements are [nF x nK x N x Nudof]
    % 
    % See also computeK, Waveguide.
    % 
    % 2022 - Daniel A. Kiefer
    % Institut Langevin, Paris, France
    % 
    if ~isvector(k), error('Wavenumbers should be a [1xN] array.'); end
    k = k(:).'; % row vector
    for i=1:length(guws) % solve for a list of waveguide objects
        guw = guws(i);
        if nargin < 3, nModes = size(guw.op.M,1); end % default
        if nModes > size(guw.op.M,1)
            warning('GEWTOOL:computeW:tooManyModes', 'More modes requested than available. Resetting nModes to the matrix size.')
            nModes = size(guw.op.M,1);
        end
        kh = k*guw.np.h0;
        M = guw.op.M; L0 = guw.op.L0; L1 = guw.op.L1; L2 = guw.op.L2;
        whn = nan(nModes, length(kh));
        u = nan(nModes, length(kh), guw.geom.Ndof);
        gdoffree = setdiff([guw.geom.gdofOfLay{:}], guw.geom.gdofDBC(:).');
        for n = 1:length(kh)
            [un, whn2] = polyeig((1i*kh(n))^2*L2 + (1i*kh(n))*L1 + L0, M); % does not work properly with eig()
            whnn = real(sqrt(whn2));
            spurious = whnn==0; whnn(spurious) = nan; un(:,spurious) = nan;
            [whnn, ind] = sort(whnn);
            un = un(:,ind);
            u(:,n,gdoffree) = un(:,1:nModes).'; % save
            whn(:, n) = whnn(1:nModes);
        end
        % save to output variable:
        dat(i).w = whn*guw.np.fh0/guw.np.h0;
        dat(i).k = kh.*ones(size(whn))/guw.np.h0;
        dat(i).u = cell(guw.geom.nLay, 1);
        for l = 1:guw.geom.nLay
            ulay = u(:,:,guw.geom.gdofOfLay{l});
            dat(i).u{l} = reshape(ulay, [size(whn), guw.geom.N(l), guw.geom.Nudof(l)]);
        end
    end
end
