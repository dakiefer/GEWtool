function dat = computeK(gews, w, nModes, opts)
    % computeK - Obtain complex wavenumbers k for specified frequencies w.
    % Solves the polynomial eigenvalue problem [(ik)^2*L2 + ik*L1 + L0(w)]*u = 0.
    %
    % Arguments:
    % - gews:    Waveguide object(s), either a scalar or vector.
    %            Describes the eigenproblem, i.e., the matrices Li.
    %            If gews is a vector, computeW solves one problem after another 
    %            and returns a vector of results "dat" of same length.
    % - w:       Angular frequencies to specify in rad/s. Vector valued.
    % - nModes:  (optional) Number of modes to compute/save (discards the highest wavenumbers).
    % - opts:    (optional) A structure of options. Possible fields are: 
    %            - 'sparse': false (default) | true. Use sparse matrices.
    %            - 'subspace': false (default) | true. Use eigs() instead of eig().
    %            - 'parallel': false (default) | true. Multi-core computation.
    %
    % Return value:
    % - dat:     A data structure containing 
    %            - w: the angular frequencies in rad/s, expanded to [nK x nF]
    %            - k: the wavenumbers in rad/m [nK x nF]
    %            - u: the displacement eigenvectors as a 
    %                 cell array describing the layers, elements are [nK x nF x N x Nudof]
    % 
    % See also computeW, Waveguide.
    % 
    % 2022-2023 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France

    if nargin < 4, opts = []; end
    if ~isfield(opts, 'sparse'),   opts.sparse = false;    end
    if ~isfield(opts, 'subspace'), opts.subspace = false;  end
    if ~isfield(opts, 'parallel')  % if a parallel pool already exists, then use it!
        parPool = gcp('nocreate'); % is empty if no parallel pool exists
        if isempty(parPool), opts.parallel = false; else, opts.parallel = true; end
    end
    if ~opts.subspace && opts.sparse
        warning('GEWTOOL:computeK:ignoringSparse',...
            'Sparse matrices are only supported in combination with the subspace solver, i.e., eigs(). Switching to subspace method. To hide this message set opts.subspace=true;');
        opts.subspace = true; % switch to eigs()
    end 
    if opts.parallel, opts.parallel = inf; else, opts.parallel = 0; end % set the number of workers
    if nargin >= 3 && ~isempty(nModes) && ~isinf(nModes) && ( ~isscalar(nModes) || nModes ~= round(nModes) )
        error('GEWTOOL:wrongArg', 'Argument "nModes" must be a scalar integer.');
    end
    
    if ~isvector(w), error('Angular frequencies should be a [Nx1] array.'); end
    w = w(:).'; % row vector
    for i = 1:length(gews) % solve for a list of waveguide objects
        gew = gews(i);
        wh = w*gew.np.h0;
        if opts.sparse
            M = sparse(gew.op.M); L0 = sparse(gew.op.L0); L1 = sparse(gew.op.L1); L2 = sparse(gew.op.L2);
        else
            M = gew.op.M; L0 = gew.op.L0; L1 = gew.op.L1; L2 = gew.op.L2;
        end
        if nargin < 3 || isempty(nModes) || isinf(nModes)
            nModes = size(gew.op.M,1);  % for now we dont distinguish between linearized and quadratic EVP
        end
        if nModes > size(gew.op.M,1)
            warning('GEWTOOL:computeK:tooManyModes', 'More modes requested than available. Resetting nModes to the matrix size.')
            nModes = size(gew.op.M,1);
        end
        kh = nan(nModes, length(wh));
        u = zeros(nModes, length(wh), gew.geom.Ndof);
        gdoffree = gew.geom.gdofFree;
        parfor (ii = 1:length(wh), opts.parallel)
            whn = wh(ii)/gew.np.fh0; % current frequency-thickness (normalized)
            [un, khn] = solveAtFreq(L2, L1, L0, M, whn, nModes, gew.geom, opts);
            u(:,ii,gdoffree) = un(:,1:nModes).'; % save
            kh(:,ii) = khn(1:nModes);
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

function [un, khn] = solveAtFreq(L2, L1, L0, M, whn, nModes, geom, opts)
    if all(L2 == 0, 'all') % is linearized as [(ik) L1 + L0 + w^2 M]*u = 0 (is this ever used?)
        if opts.subspace
            [un, khn] = eigs(L0 + whn^2*M, -1i*L1, nModes, "smallestabs");
            khn = diag(khn); % eigs returns a matrix
        else
            [un, khn] = eig(L0 + whn^2*M, -1i*L1, 'vector');
        end
%         [un, ikhi] = polyeig(L0 + whn^2*M, L1); 
    elseif isempty(L1) % is linearized as [(ik)^2 L2 + L0 + w^2 M]*u = 0
        if opts.subspace
            [un, khn2] = eigs(L0 + whn^2*M, L2, nModes, "smallestabs");
            khn = sqrt(diag(khn2)); % eigs returns a matrix
        else
            [un, khn2] = eig(L0 + whn^2*M, L2, 'vector');
            khn = sqrt(khn2);
        end
        dofy = geom.gdofRedY;
        un(dofy,:) = -1i*un(dofy,:)./khn.'; % eig.vec. was [ux, 1i*k*uy]
    else % quadratic EVP: [(ik)^2 L2 + (ik) L1 + L0 + w^2 M]*u = 0
        if opts.subspace
            [un, ikhn] = GEWpolyeig(L0 + whn^2*M, L1, L2, nModes, opts);
            khn = -1i*ikhn; % extract kh
        else
            [un, ikhn] = polyeig(L0 + whn^2*M, L1, L2); % compute 1i*kh: double as fast than computing kh
            khn = -1i*ikhn; % extract kh
        end
    end
    spurious = isinf(khn) | isnan(khn); khn(spurious) = []; un(:,spurious) = [];
    [khn, ind] = sort(khn);
    un = un(:,ind); % sort
end
