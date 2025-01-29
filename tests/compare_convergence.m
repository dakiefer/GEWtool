%% Compare the convergence of different discretization methods
% Compares accuracy vs matrix size of the following methods: 
% - FEM: linear finite elements
% - SCM: spectral collocation method (Chebyshev)
% - SEM (GEWtool): spectral element on GLL nodes and Lagrange polynomials 
%
% The test is randomly done for the 4th mode (S1 mode) at 5.6 rad/m of a 1-m thick
% aluminum plate (approx. 5 kHz).
% 
% The error is quantified by comparing against the Rayleigh-Lamb root computed
% with numerical precision with Matlab's fzero() function.
% 
% 2025 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France

mat = Material('aluminum'); % load material parameters from database
indMode = 4;    % we look at the 4th mode
w0 = 2*pi*5000; % for reference and mode selection
k0 = 5.6;       % for reference and mode selection (also: 3.6180, 3.8674, 5.59596, 7.8648)
Ns = 4:40;       % list of discretization orders
solver = @(L2, L1, L0, M) eigs(k0^2*L2 - 1i*k0*L1 - L0, M, 8, 'smallestabs');
tolw = 1e-10; % we consider solutions to be accurate with this tolerance
 
% % compute Rayleigh-Lamb root:
rayLambAtK0 = @(wh) rayleighLambS(mat, wh, k0).*rayleighLambA(mat, wh, k0);
% options = optimset('Display','iter'); % uncomment to show convergence
options.TolX = eps;
wRayLamb = fzero(rayLambAtK0, w0, options); % only for real arguments
residuumAtRLRoot = abs(rayLambAtK0(wRayLamb));
assert(residuumAtRLRoot < 10e-14); % make sure that Rayleigh-Lamb root is converged before starting tests

% % % SEM with GEWtool
generateMats = @(N) matrices_SEM(mat, N);
[wsSEM, nsSEM] = computeFrequencies(generateMats, solver, Ns, indMode, w0); 
[errSEM, nConvSEM] = errorAnalysis(wsSEM, nsSEM, wRayLamb, tolw); 

% % FEM (linear) with GEWtool
generateMats = @(N) matrices_FEM(mat, N);
[wsFEM, nsFEM] = computeFrequencies(generateMats, solver, Ns, indMode, w0); 
[errFEM, nConvFEM] = errorAnalysis(wsFEM, nsFEM, wRayLamb, tolw); 

% % % SCM
generateMats = @(N) matrices_SCM(mat, N);
[wsSCM, nsSCM] = computeFrequencies(generateMats, solver, Ns, indMode, w0); 
[errSCM, nConvSCM] = errorAnalysis(wsSCM, nsSCM, wRayLamb, tolw); 

%% plot
figure(1), clf; hold on, 
plot(nsFEM, errFEM, 'o-', 'Color', 0.6*[1 1 1], 'DisplayName', 'FEM');
plot(nsSCM, errSCM, 'o-', 'DisplayName', 'SCM', 'SeriesIndex',2);
plot(nsSEM, errSEM, 'o-', 'DisplayName', 'SEM', 'SeriesIndex',1);
ax=gca; ax.YScale='log';
c = mat.c; rho = mat.rho;
% yline(tolw)
% plot(nConvSCM, errSCM(nsSCM == nConvSCM), 'ro');
xlabel('matrix size $n$'), ylabel('relative error $(\omega - \omega_0)/\omega_0$')
title(sprintf('S1 mode at $k h$ = %g rad', k0))
legend('Location','southwest'); 

disp('SCM converges with')
disp(nConvSCM)
disp('SEM converges with')
disp(nConvSEM)



% % HELPER FUNCTIONS:
% % to assembling the matrices.
function [ws, ns] = computeFrequencies(generateMats, solver, Ns, indMode, w0)
    ws = nan(numel(Ns),1); % allocate
    ns = nan(numel(Ns),1); % allocate for resulting matrix size (depends on considered polarization)
    for i=1:numel(Ns)
        [L2, L1, L0, M, f0] = generateMats(Ns(i));
        % time(i) = timeit(@()generateMats(N(i))); % test assembling time
        [w2] = solver(L2, L1, L0, M);  % use Cholesky: positive definite B
        w = sort(sqrt(w2))*f0;
        % [~, indSel] = min(abs(w - w0));
        ws(i)=w(indMode);
        ns(i)=size(L0,1);
    end
end

function [errRel, nConv] = errorAnalysis(ws, ns, wRayLamb, tolw)
    errRel = abs(ws-wRayLamb)/wRayLamb;
    converged = abs(errRel) < tolw;        % indices of converged values
    indConv = find( diff(converged) ) + 1;  % indices where the drops occur
    % assert( length(indConv) == 1 );     % there exists exactly one index where the drop occurs
    nConv = ns(indConv); 
end

function [L2, L1, L0, M, f0] = matrices_SEM(mat, N)
    % matrices_SEM - one element with N nodes (of order N-1)
    plate = Plate(mat, 1, N);
    plate.np.h0 = 1; 
    gew = plate.Lamb; 
    L2 = gew.op.L2; L1 = gew.op.L1; L0 = gew.op.L0; M = gew.op.M; 
    f0 = gew.np.fh0;
end

function [L2, L1, L0, M, f0] = matrices_FEM(mat, N)
    % matrices_FEM - Finite Elements of order Nnodes-1
    h = 1;
    Nnodes = 2; 
    Nelem = N-1; 
    mats = repmat(mat, 1, Nelem);
    plate = Plate(mats, h/Nelem, Nnodes);
    plate.np.h0 = 1; % important here! otherwise h0 changes with Nelem
    gew = plate.Lamb; 
    L2 = gew.op.L2; L1 = gew.op.L1; L0 = gew.op.L0; M = gew.op.M; 
    f0 = gew.np.fh0;
end

function [L2, L1, L0, M, f0] = matrices_SCM(mat, N)
    c = mat.c; rho = mat.rho;
    c0 = c(1,2,1,2); rho0 = rho; f0 = sqrt(c0/rho0);  % normalization parameters
    rhon = rho/rho0; cn = c/c0; % normalize
    
    % relevant material matrices: 
    udof = 1:2; % Lamb [1:2]; SH [3]; coupled [1:3];
    cxx = squeeze(cn(1,udof,udof,1));
    cxy = squeeze(cn(1,udof,udof,2)); 
    cyx = squeeze(cn(2,udof,udof,1));
    cyy = squeeze(cn(2,udof,udof,2));
    I = eye(size(cxx)); 
    
    % discretization 
    [~, Dy_dash] = chebdif(N, 2); % create differentiation matrices
    D1 = -2*Dy_dash(:,:,1); % differentiation on unit domain
    D2 = 4*Dy_dash(:,:,2);  % second order derivative
    Id = eye(size(D1));     % identity matrix for discretization
    
    % define wave operators:
    L2 = kron(cxx, Id); L1 = kron(cxy + cyx, D1); L0 = kron(cyy, D2); 
    M = kron(rhon*I, Id);
    B1 = kron(cyx, Id([1, N], :)); B0 = kron(cyy, D1([1, N], :));
    
    % incorporate BCs:
    dofBC = [(0:length(udof)-1)*N+1; (1:length(udof))*N]; % [1, N, N+1, 2*N, 2*N+1, 3*N];
    L2(dofBC, :) = 0; L1(dofBC, :) = B1; L0(dofBC, :) = B0; M(dofBC, :) = 0;
end

function [L2, L1, L0, M, f0] = matrices_SCMrect(mat, N)
    c = mat.c; rho = mat.rho;
    c0 = c(1,2,1,2); rho0 = rho; f0 = sqrt(c0/rho0); % normalization parameters
    rhon = rho/rho0; cn = c/c0;
    
    % relevant stiffness tensors: 
    cxx = squeeze(cn(1,1:2,1:2,1));
    cxy = squeeze(cn(1,1:2,1:2,2));
    cyx = squeeze(cn(2,1:2,1:2,1));
    cyy = squeeze(cn(2,1:2,1:2,2));
    Rho = rhon*eye(size(cxx));
    
    % % discretization 
%     [~, Dy_dash] = chebdif(N, 2);
%     D1 = -2*Dy_dash(:,:,1); % differentiation on [-1/2, 1/2]
%     D2 = (2)^2*Dy_dash(:,:,2); % 2nd order differentiation on [-1/2, 1/2]
    D1 = diffmat(N, 1, [-0.5 0.5], 'chebkind2'); % first order on domain [0 1]
    D2 = diffmat(N, 2, [-0.5 0.5], 'chebkind2'); % second order on domain [0 1]
    Id = eye(size(D1));  % identity matrix for discretization
    [xf, ccw, baryw] = chebpts(N, [-0.5 0.5], 2); % second-kind points (includes boundaries)
    [xd] = chebpts(N-2, [-0.5 0.5], 1); % first-kind points (without boundaries)
    P = barymat(xd, xf, baryw); % resampling matrix (from values on xf to values on xd)
    % P = eye(size(Id)); % disable resampling
    
    % % Lamb wave problem
    L2 = kron(cxx, P*Id); L1 = kron(cxy + cyx, P*D1); L0 = kron(cyy, P*D2); 
    M = kron(Rho,P*Id);
    B1 = kron(cyx, Id([1, N], :)); B0 = kron(cyy, D1([1, N], :)); % BCs
    L2 = [zeros(size(B0)); L2]; L1 = [B1; L1]; L0 = [B0; L0]; M = [zeros(size(B0)); M]; % incorporate BCs
end

