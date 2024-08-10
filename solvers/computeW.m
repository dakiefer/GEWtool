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
    %            - 'eigenvecs': true (default) | false. Whether to compute
    %               eigenvectors.
    %            - 'standardEVP': true | false. Whether to convert the
    %               generalized eigenvalue problem to a standard one. Default: true
    %               if M is diagonal, false otherwise.
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
    % 2022-2024 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France

    if nargin < 4, opts = []; end
    if nargin < 3, nModes = []; end
    
    if ~isvector(k), error('Wavenumbers should be a [1xN] array.'); end
    k = k(:); % column vector
    for i=1:length(gews) % solve for a list of waveguide objects
        gew = gews(i);
        opti = parseSolverOpts(opts,gew.op,nModes); % opti might be modified in the iteration
        if nargin < 3 || isempty(nModes) || isinf(nModes) % default
            nModes = size(gew.op.M,1); 
        end
        if nModes > size(gew.op.M,1)
            warning('GEWTOOL:computeW:tooManyModes', 'More modes requested than available. Resetting nModes to the matrix size.')
            nModes = size(gew.op.M,1);
        end
        if isfield(opti,'target') & isnumeric(opti.target)
            opti.target = (opti.target*gew.np.h0/gew.np.fh0)^2;
        else 
            opti.target = 'smallestabs'; % parfor throws error if target undefined
        end
        if opti.sparse
            M = sparse(gew.op.M); L0 = sparse(gew.op.L0); L1 = sparse(gew.op.L1); L2 = sparse(gew.op.L2);
        else
            M = gew.op.M; L0 = gew.op.L0; L1 = gew.op.L1; L2 = gew.op.L2;
        end
        if opti.standardEVP
            [H2, H1, H0, T] = transformToStandardEVP(L2, L1, L0, M);
            solveAtK = @(kh) solveEVP(kh^2*H2 - 1i*kh*H1 - H0, nModes, opti);
            opti.T = T; 
        else
            solveAtK = @(kh) solveGEP(kh^2*L2 - 1i*kh*L1 - L0, M, nModes, opti);
        end
        gdoffree = gew.geom.gdofFree;
        kh = k*gew.np.h0;
        whn = nan(length(kh), nModes); % allocate
        if opti.eigenvecs
            u = zeros(length(kh), nModes, gew.geom.Ndof); % allocate
            parfor (j = 1:length(kh), opti.parallel)
                [lbd, eVec] = solveAtK(kh(j));
                [whnj, uj] = retrieveWu(lbd, eVec, nModes, opti); 
                whn(j,:) = whnj;
                u(j,:,gdoffree) = uj.'; % save
            end
        else
            parfor (j = 1:length(kh), opti.parallel)
                lbd = solveAtK(kh(j));
                whn(j,:) = retrieveW(lbd, nModes);
            end
        end
        % save to output variable:
        if ~gew.isDissipative
            whn = real(whn); % remove imaginary part due to numerical inaccuracy
        end
        dat(i).w = whn*gew.np.fh0/gew.np.h0;
        dat(i).k = kh.*ones(size(whn))/gew.np.h0;
        if opti.eigenvecs
            dat(i).u = cell(gew.geom.nLay, 1);
            for l = 1:gew.geom.nLay
                ulay = u(:,:,gew.geom.gdofOfLay{l});
                dat(i).u{l} = reshape(ulay, [size(whn), gew.geom.N(l), gew.geom.Nudof(l)]);
            end
        end
    end
end

function whn = retrieveW(lbd, nModes)
    whn = sort(sqrt(lbd));
    whn = whn(1:nModes);
end

function [whn, un] = retrieveWu(lbd, eVec, nModes, opti)
    [whn, ind] = sort(sqrt(lbd));
    whn = whn(1:nModes);
    eVec = eVec(:,ind); eVec = eVec(:,1:nModes);
    if opti.standardEVP % transform eigenvectors back to displacements
        eVec = opti.T*eVec;
    end
    un = eVec;
end

function [H2, H1, H0, T] = transformToStandardEVP(L2, L1, L0, M)
    % transformToStandardEVP - transform to a standard eigenvalue problem
    % (-ik^2*L2 - ik*L1 - L0)*u = lbd*M*u   ->   (-ik^2*H2 - ik*H1 - H0)*y = lbd*I*y
    % Using T = M^(-1/2), we transform into Hi = T*Li*T. The original eigenvectors 
    % can be recovered from u = T*y. 
    %
    % Arguments:
    % - L2, L1, L0, M: matrices of the original eigenvalue problem   
    %
    % Return values:
    % - H2, H1, H0:   Matrices of the standard eigenvalue problem (EVP)
    % - T:            Matrix that transforms eigenvectors y to u, i.e., T*y = u
    % 
    % 2024 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France
    %        Malte Röntgen, LAUM, Le Mans Université, France
    if ~issparse(M) && isinf(cond(M)) || issparse(M) && isinf(condest(M))
        error('GEWTOOL:transformToStandardEVP', 'The mass matrix must be regular in order to transform the generalized eigenvalue problem to a standard one. Set "opts.standardEVP = false" and pass "opts" to computeW() in order to deactivate the transformation. See "help computeW" for more information.');
    elseif ~issparse(M) && cond(M) > 1e4 || issparse(M) && condest(M) > 1e4
        warning('GEWTOOL:transformToStandardEVP', 'The mass matrix is badly conditioned. Consider setting "opts.standardEVP = false" and pass "opts" to computeW() in order to deactivate the transformation. See "help computeW" for more information.')
    end
    if isdiag(M)
        if issparse(M)
            T = sparse(diag(1./sqrt(diag(M))));
        else
            T = diag(1./sqrt(diag(M))); % usually the case in GEWtool
        end
    else
        T = M^(-1/2); 
    end
    H2 = T*L2*T; H1 = T*L1*T; H0 = T*L0*T;
    % recover Hermitean symmetry:
    if ishermitian(L2) && ishermitian(1i*L1) && ishermitian(L0) && ishermitian(T)
        H2 = (H2 + H2')/2;
        H1 = (H1 - H1')/2; % skew-hermitian
        H0 = (H0 + H0')/2;
    end
end
