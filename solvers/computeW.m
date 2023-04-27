function dat = computeW(gews, k, nModes, opts)
    % computeW - Obtain frequency w for specified wavenumbers k.
    % Solves the eigenvalue problem described by -L(k)*u = w^2*M*u.
    %
    % Arguments:
    % - gews:    Waveguide object(s), either a scalar or vector.
    %            Describes the eigenproblem, i.e., the matrices L and M.
    %            If gews is a vector, computeW solves one problem after another 
    %            and returns a vector of results "dat" of same length.
    % - k:       Wavenumbers to specify in rad/m. Vector valued.
    % - nModes:  (optional) Number of modes to compute/save (discards the highest frequencies).
    % - opts:    (optional) A structure of options. Possible fields are: 
    %            - 'sparse': false (default) | true. Use sparse matrices.
    %            - 'subspace': false (default) | true. Use eigs() instead of eig().
    %            - 'parallel': false (default) | true. Multi-core computation.
    %
    % Return value:
    % - dat:     A data structure containing 
    %            - w: the angular frequencies in rad/s [nK x nF]
    %            - k: the wavenumbers in rad/m, expanded to [nK x nF]
    %            - u: the displacement eigenvectors as a 
    %                 cell array describing the layers, elements are [nK x nF x N x Nudof]
    % 
    % See also computeK, Waveguide.
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
        warning('GEWTOOL:computeW:ignoringSparse',...
            'Sparse matrices are only supported in combination with the subspace solver, i.e., eigs(). Switching to subspace method. To hide this message set opts.subspace=true;');
        opts.subspace = true; % switch to eigs()
    end 
    if opts.parallel, opts.parallel = inf; else, opts.parallel = 0; end % set the number of workers
    if nargin >= 3 && ~isempty(nModes) && ~isinf(nModes) && ( ~isscalar(nModes) || nModes ~= round(nModes) )
        error('GEWTOOL:wrongArg', 'Argument "nModes" must be a scalar integer.');
    end
    
    if ~isvector(k), error('Wavenumbers should be a [1xN] array.'); end
    k = k(:); % column vector
    for i=1:length(gews) % solve for a list of waveguide objects
        gew = gews(i);
        if nargin < 3 || isempty(nModes) || isinf(nModes) % default
            nModes = size(gew.op.M,1); 
        end
        if nModes > size(gew.op.M,1)
            warning('GEWTOOL:computeW:tooManyModes', 'More modes requested than available. Resetting nModes to the matrix size.')
            nModes = size(gew.op.M,1);
        end
        kh = k*gew.np.h0;
        if opts.sparse
            M = sparse(gew.op.M); L0 = sparse(gew.op.L0); L1 = sparse(gew.op.L1); L2 = sparse(gew.op.L2);
        else
            M = gew.op.M; L0 = gew.op.L0; L1 = gew.op.L1; L2 = gew.op.L2;
        end
        whn = nan(length(kh), nModes);
        u = zeros(length(kh), nModes, gew.geom.Ndof);
        gdoffree = setdiff([gew.geom.gdofOfLay{:}], gew.geom.gdofDBC(:).');
        useSubspace = opts.subspace; % extracting option avoids Matlab warning due to parfor loop
        parfor (n = 1:length(kh), opts.parallel)
            if useSubspace
                [un, whn2] = eigs(-(1i*kh(n))^2*L2 - (1i*kh(n))*L1 - L0, M, nModes, "smallestabs");
                whn2 = diag(whn2); % eigs returns a matrix
            else
                [un, whn2] = eig(-(1i*kh(n))^2*L2 - (1i*kh(n))*L1 - L0, M, 'chol',...
                'vector'); % Choleski guarantees real eigenvalues, M needs to be positive definite, faster than qz
            end
            [whnn, ind] = sort(sqrt(whn2));
            un = un(:,ind);
            whn(n,:) = real(whnn(1:nModes)); % save
            u(n,:,gdoffree) = un(:,1:nModes).'; % save
        end
        % save to output variable:
        dat(i).w = whn*gew.np.fh0/gew.np.h0;
        dat(i).k = kh.*ones(size(whn))/gew.np.h0;
        dat(i).u = cell(gew.geom.nLay, 1);
        for l = 1:gew.geom.nLay
            ulay = u(:,:,gew.geom.gdofOfLay{l});
            dat(i).u{l} = reshape(ulay, [size(whn), gew.geom.N(l), gew.geom.Nudof(l)]);
        end
    end
end
