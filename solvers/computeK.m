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
    %            - 'eigenvecs': true (default) | false. Whether to compute
    %               eigenvectors.
    %            - 'standardEVP': true | false. Whether to convert the
    %               generalized eigenvalue problem to a standard one. Default: true
    %               if M is diagonal, false otherwise.
    %            - 'sparse': false (default) | true. Use sparse matrices.
    %            - 'subspace': false (default) | true. Use eigs() instead of eig().
    %            - 'parallel': false (default) | true. Multi-core computation.
    %            - 'show': print the used options when computing (for debugging)
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
    % 2022-2024 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France

    if nargin < 4, opts = []; end
    if nargin < 3, nModes = []; end
    
    if ~isvector(w), error('Angular frequencies should be a [Nx1] array.'); end
    w = w(:).'; % row vector
    dat = repmat(GEWdat(gews(1),[],w,[]),1,length(gews));
    for i = 1:length(gews) % solve for a list of waveguide objects
        gew = gews(i);
        [opti, nModes] = parseSolverOpts(opts, gew.op, nModes); % opti will be modified in the iteration
        if isfield(opti,'target') & isnumeric(opti.target)
            opti.target = opti.target*gews.np.h0; % normalize to match gew.op
        end
        if opti.sparse
            M = sparse(gew.op.M); L0 = sparse(gew.op.L0); L1 = sparse(gew.op.L1); L2 = sparse(gew.op.L2);
        else
            M = gew.op.M; L0 = gew.op.L0; L1 = gew.op.L1; L2 = gew.op.L2;
        end
        if isempty(L1)
            if isfield(opti,'target') & isnumeric(opti.target)
                opti.target = opti.target^2;
            end
            A = L0; B = L2; AA = M; 
            opti.linearization = 'k2';
        else
            if isfield(opti,'target') & isnumeric(opti.target)
                opti.target = 1i*opti.target;
            end
            [A, B, AA] = linearizePolyEig(L2, L1, L0, M);
            opti.linearization = 'companion';
        end
        n = size(L0,1);
        if opti.standardEVP
            [H, HH, T] = transformToStandardEVP(A, B, AA);
            solveAtW = @(wh) solveEVP(H + wh^2*HH, nModes, opti);
            switch opti.linearization 
                case 'companion'
                    secondBlockInd = (n+1):2*n;
                    opti.T = T(secondBlockInd,:); % eigenvector is [ik*u, u]
                case 'k2'
                    opti.T = T;
            end
        else
            solveAtW = @(wh) solveGEP(A + wh^2*AA, B, nModes, opti);
        end
        geom = gew.geom; 
        gdoffree = geom.gdofFree;
        whn = w*gew.np.h0/gew.np.fh0;
        kh = nan(nModes, length(whn));
        if opti.eigenvecs
            u = zeros(nModes, length(whn), length(gdoffree)); % allocate
            parfor (j = 1:length(whn), opti.parallel)
                [lbd, eVec] = solveAtW(whn(j));
                [khj, uj] = retrieveKu(lbd, eVec, nModes, opti, geom);
                kh(:,j) = khj; 
                u(:,j,:) = uj.';
            end
        else
            parfor (j = 1:length(whn), opti.parallel)
                lbd = solveAtW(whn(j));
                kh(:,j) = retrieveK(lbd, nModes, opti);
            end
            u = [];
        end
        k = kh/gew.np.h0;
        dat(i) = GEWdat(gew,k,w,u); % save in an object of class 'GEWdat' 
    end
end

function khn = retrieveK(lbd, nModes, opti)
    sortAccuracy = 1e6;
    switch opti.linearization
        case 'companion'
            khn = -1i*lbd;
        case 'k2'
            khn = sqrt(lbd);
    end
    khnRounded = round(khn*sortAccuracy)/sortAccuracy; % sort on digits with sufficient presition only
    [~, ind] = sort(khnRounded,'ComparisonMethod','abs'); % sort by real part, 
    khn = khn(ind); 
    khn = khn(1:nModes); % sort and save kh(:,n) = 
end

function [khn, un] = retrieveKu(lbd, eVec, nModes, opti, geom)
    sortAccuracy = 1e6;
    switch opti.linearization
        case 'companion'
            khn = -1i*lbd;
        case 'k2'
            khn = sqrt(lbd);
            dofz = geom.gdofRedZ;
            eVec(dofz,:) = -1i*eVec(dofz,:)./khn.'; % eig.vec. was [ux, 1i*k*uy]
    end
    khnRounded = round(khn*sortAccuracy)/sortAccuracy; % sort on digits with sufficient presition only
    [~, ind] = sort(khnRounded,'ComparisonMethod','abs'); % sort by real part, 
    khn = khn(ind); khn = khn(1:nModes); % sort and crop
    eVec = eVec(:,ind); eVec = eVec(:,1:nModes); % sort and crop
    if     strcmp(opti.linearization,'companion') &&  opti.standardEVP 
        eVec = opti.T*eVec; % transform eigenvectors back 
    elseif strcmp(opti.linearization,'companion') && ~opti.standardEVP
        n = length(eVec)/2;
        eVec = eVec(n+1:2*n,:); % eigenvectors are [ik*u, u] 
    end
    un = eVec;
end

function [A, B, AA] = linearizePolyEig(L2, L1, L0, M)
    % linearizePolyEig - companion linearization of polynomial eigenvalue problem.
    % (ik^2*L2 + ik*L1 + L0 + w^2*M)*u = 0    ->    A*x = ik*B*x
    % 
    % See also polyeig.
    % 
    % 2024 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France
    %        Malte Röntgen, LAUM, Le Mans Université, France

    % % companion matrix linearization
    p = 2; % polynomial order (only two for now)
    n = size(L0,1); 
    firstBlockInd = 1:n;
    secondBlockInd = (n+1):(2*n);
    
    % allocate matrices:
    if issparse(L0)
        A = sparse(n*p, n*p);
        B = sparse(n*p, n*p); % sparse zeros-matrix
        AA = sparse(n*p, n*p);
    else
        A = zeros(n*p);
        B = zeros(n*p);
        AA = zeros(n*p);
    end
    
    % First companion linearization (A + w*AA - lbd*B):
    % construct A:
    A(firstBlockInd,secondBlockInd)  = eye(n);
    A(secondBlockInd,secondBlockInd) = -L1; 
    A(secondBlockInd,firstBlockInd)  = -L0;
    % construct B: (NOTE this could be avoided when converting to standard EVP since
    % only B^(-1/2) is needed)
    B(firstBlockInd,firstBlockInd) = eye(n); 
    B(secondBlockInd,secondBlockInd) = L2;
    % construct AA: 
    AA(secondBlockInd,firstBlockInd)  = -M;
end

function [H, HH, T] = transformToStandardEVP(A, B, AA)
    % transformToStandardEVP - transform to a standard eigenvalue problem
    % (A + w^2*AA)*u = ik*B*u   ->   (H + w^2*HH)*y = ik*I*y
    % Using T = B^(-1/2), we transform into H = T*A*T and HH = T*AA*T. The
    % original eigenvectors can be recovered from u = T*y.
    %
    % Arguments:
    % - A, B, AA: matrices of the original eigenvalue problem
    %
    % Return values:
    % - H, HH:    Matrices of the standard eigenvalue problem (EVP)
    % - T:        Matrix that transforms eigenvectors y to u, i.e., T*y = u
    %
    % 2024 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France
    %        Malte Röntgen, LAUM, Le Mans Université, France
    if isdiag(B)
        if issparse(B)
            T = sparse(diag(1./sqrt(diag(B))));
        else
            T = diag(1./sqrt(diag(B))); % usually the case in GEWtool
        end
    else
        T = B^(-1/2); 
    end
    H = T*A*T; HH = T*AA*T; 
    % recover Hermitean symmetry:
    if ishermitian(A) && ishermitian(AA) && ishermitian(T)
        H = (H + H')/2;
        HH = (HH + HH')/2;
    end
end
