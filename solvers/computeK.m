function dat = computeK(gews, w, nModes, opts)
% computeK - Obtain complex wavenumbers k for specified frequencies w.
% Solves the polynomial eigenvalue problem [(ik)^2*L2 + ik*L1 + L0(w)]*u = 0.
%
% Arguments:
% - gews:    Waveguide object(s), either a scalar or vector.
%            Describes the eigenproblem(s), i.e., the matrices Li.
%            If gews is a vector, computeW solves one problem after another 
%            and returns a vector of results "dat" of the same length.
% - w:       Angular frequencies to specify in rad/s. Vector valued.
% - nModes:  (optional) Number of modes to compute/save (discards highest).
% - opts:    (optional) A structure of options. Possible fields are: 
%            - 'eigenvecs': true (default) | false. Whether to compute eigenvectors.
%               Turn off for speedup. 
%            - 'standardEVP': true | false. Whether to convert the
%               generalized eigenvalue problem to a standard one. Default: true
%               if M is diagonal, false otherwise.
%            - 'subspace': false | true. Use eigs() instead of eig() for speedup.
%               Defaults to true when size(op.M,1) > 60. You should also provide 
%               nModes to computeK.
%            - 'sparse': false | true. Use sparse matrices. Default == 'subspace'.
%            - 'parallel': false | true. Multi-core computation. Defaults to true 
%               when a parallel pool is running.
%            - 'trace': true (default) | false. Reorder modes so that they are properly
%               sorted (grouped). Turn off for consistency across multiple computations. 
%            - 'show': print the used options when computing (for debugging)
%
% Return value:
% - dat:     Object(s) of class 'GEWdat' that stores 
%            - w: the angular frequencies in rad/s, expanded to [nK x nW]
%            - k: the wavenumbers in rad/m [nK x nW]
%            - Psi: the eigenvectors [nK x nW x size(gews(i).op.M,2)]
%            - u: the (displacement) eigenvectors expanded to a cell array where
%            each entry is the eigenvector on a layer of size [nK x nW x Nlay x Nudof]
% 
% See also computeW, Waveguide.
% 
% 2022-2026 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France

    if nargin < 4, opts = struct(); end
    if nargin < 3, nModes = []; end
    opts.compute = "k";
    
    if ~isvector(w), error('Angular frequencies should be a [Nx1] array.'); end
    w = w(:).'; % row vector
    dat = repmat(GEWdat(gews(1),[],w,[]),1,length(gews));
    for i = 1:length(gews) % solve for a list of waveguide objects
        gew = gews(i);
        opts.solver = @computeK;
        [opti, nModes] = parseSolverOpts(opts, gew.op, nModes); % opti will be modified in the iteration
        if isfield(opti,'target') & isnumeric(opti.target)
            opti.target = opti.target*gew.np.h0; % normalize to match gew.op
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
            if opti.eigenvecs % T not needed if only eigenvalues are computed
                switch opti.linearization
                    case 'companion'
                        secondBlockInd = (n+1):2*n;
                        opti.T = T(secondBlockInd,:); % eigenvector is [ik*u, u]
                    case 'k2'
                        opti.T = T;
                end
            end
        else
            solveAtW = @(wh) solveGEP(A + wh^2*AA, B, nModes, opti);
        end
        whn = w*gew.np.h0/gew.np.fh0;
        kh = nan(nModes, length(whn));
        if opti.eigenvecs
            u = zeros(nModes, length(whn), size(gew.op.L0,1)); % allocate
            geom = gew.geom; % extract before "parfor" (overhead due to broadcast) 
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

        if opti.trace % trace modes: re-arrange such that dat.k
            ind = []; % initialize
            try
                [kh, ind] = reorderByProximity(kh); % matches the modes such that the wavenumbers change as little as possible with frequency
            catch exception
                warning('I was not able to re-order modes by proximity (tracing). Modes will be unordered.')
            end
            if ~isempty(u) && ~isempty(ind) % order eigenvectors if mode matching was successful
                for n = 1:size(kh,2), u(:,n,:) = u(ind(:,n),n,:); end % same ordering 
            end
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
    khn = khn(ind(1:nModes)); % sort and save kh(:,n)
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
    khn = khn(ind(1:nModes));
    eVec = eVec(:,ind(1:nModes)); % sort and crop
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
    % (ik^2*L2 + ik*L1 + L0 + w^2*M)*u = 0    ->    [A + w^2*AA]*x = ik*B*x
    % where the new vectors are    x = [  u   ]
    %                                  [ ik*u ]
    % and the matrices are given by
    % 
    % A = [  0    I  ]      ,     AA = [  0   0 ]
    %     [ -L0  -L1 ]                 [ -M   0 ]
    % 
    % B = [ I   0  ]   (is diagonal and regular -> important for eigs, etc.)
    %     [ 0   L2 ]
    % 
    % See also polyeig.
    % 
    % 2024 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France
    %        Malte Röntgen, LAUM, Le Mans Université, France

    if isempty(L2) % already linear
        A = -L0; AA = -M; B = L1; 
        return; 
    end

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
