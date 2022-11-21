%% plot the convergence with respect to the increase in degrees of freedom (dof)
% The convergence is tested against the Rayleigh-Lamb root, i.e., valid for Lamb waves
% in an isotropic plate. The provided approximate solution (w0, k0) is for an
% aluminum plate.

% % input: You need to create matrices Li, M first!
matrices = @(c, rho, N) matrices_GEWtool(c, rho, N); % select one of the function defined below
addpath("~/Projekte/GEWtool/tests/") % for rayleighLambS() and rayleighLambA()
mat = Material('aluminum'); % load material parameters from database
c = mat.c; rho = mat.rho;
h = 1;          % thickness 
w0 = 2*pi*5000; % for reference and mode selection
k0 = 5.6;       % for reference and mode selection (also: 3.6180, 3.8674, 5.59596, 7.8648)
N = 4:30;       % list of discretization orders
 
% % compute Rayleigh-Lamb root:
rayLambAtK0 = @(wh) rayleighLambS(mat, wh, k0).*rayleighLambA(mat, wh, k0);
options = optimset('Display','iter');
options.TolX = eps;
wRayLamb = fzero(rayLambAtK0, w0, options); % only for real arguments
residuumAtRLRoot = abs(rayLambAtK0(wRayLamb))

%% computing the frequency w
ws = nan(numel(N),1); % allocate
dofs = nan(numel(N),1); % allocate for resulting matrix size (depends on considered polarization)
time = nan(numel(N),1); % allocate
for i=1:numel(N)
    [L2, L1, L0, M, fh0] = matrices(c, rho, N(i));
    dofs(i)=size(L0,1);
    time(i) = timeit(@()matrices(c, rho, N(i))); % test assembling time
    [wh2] = eig(-(1i*k0*h)^2*L2 - (1i*k0*h)*L1 - L0, M, 'qz');  % use Cholesky: positive definite B
    w = sqrt(wh2)*fh0/h;
    [~, indSel] = min(abs(w - w0));
    ws(i)=w(indSel);
end

errRel = abs(ws-wRayLamb)/wRayLamb;
errWNmax = errRel(end)

figure(1), hold on, 
bar(dofs, time); 
title('assembling time for matrices'); 
xlabel('matrix size (= 2N)'), ylabel('time in s')

figure(2), hold on, plot(dofs,abs(errRel),'o--');
ax=gca; ax.YScale='log';
xlabel('matrix size (= 2N)'), ylabel('rel. error')
title('Error w.r.t Rayleigh-Lamb root')
% legend({'SCM', 'SEM', 'GEWtool'})


%% assembling the matrices:
function [L2, L1, L0, M, fh0] = matrices_SCM(c, rho, N)
    c0 = c(1,2,1,2); rho0 = rho; fh0 = sqrt(c0/rho0);  % normalization parameters
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

function [L2, L1, L0, M, fh0] = matrices_SEM(c, rho, N)
    c0 = c(1,2,1,2); rho0 = rho; fh0 = sqrt(c0/rho0);  % normalization parameters
    rhon = rho/rho0; cn = c/c0; % normalize
    
    % relevant material matrices: 
    udof = 1:2; % Lamb and/or SH
    cxx = squeeze(cn(1,udof,udof,1));
    cxy = squeeze(cn(1,udof,udof,2));
    cyx = squeeze(cn(2,udof,udof,1));
    cyy = squeeze(cn(2,udof,udof,2));
    I = eye(size(cxx)); 
    
    % % discretize: 
    % % using Lagrange polynomials on GLL points:
    [yi] = lobpts(N, [-1, 1]); % does only work on dom = [-1 1]!!!!
    yi = yi/2; % scale to [-1/2, 1/2]
    P = chebfun.lagrange(yi);
    Pd = diff(P);
    
    % % "element" matrices for one displacement component:
    pp = elemPP(P);         % ∫ Pi*Pj dy
    pd = elemPPd(P, Pd);    % ∫ Pi*Pj' dy
    dd = elemPdPd(Pd);      % ∫ Pi'*Pj' dy
    % % assemble for the displacement components:
    M  = kron(rhon*I,pp);
    L2 = kron(cxx, pp);
    L1 = kron(cxy, pd) - kron(cyx, pd.');
    L0 = kron(cyy, dd);
end

function [L2, L1, L0, M, fh0] = matrices_GEWtool(c, rho, N)
    mat.c = c; mat.rho = rho; h = 1;
    plate = Plate(mat, h, N);
    gew = plate.Lamb; 
    L2 = gew.op.L2; L1 = gew.op.L1; L0 = gew.op.L0; M = gew.op.M; 
    fh0 = gew.np.fh0;
end

%% element matrices:
function pp = elemPP(P) 
    % elemPP - integral ∫P*Pdy of basis functions P (element mass)
    N = size(P,2);
    pp = zeros(N);
    for i = 1:N
        for j = i:N
            pp(i,j) = sum(P(:,i)*P(:,j));
            pp(j,i) = pp(i,j);
        end
    end
end

function pd = elemPPd(P, Pd) 
    % elemPPd - integral ∫P*P'dy of basis functions P (element stiffness and flux)
    N = size(P,2);
    pd = zeros(N);
    for i = 1:N
        for j = 1:N
            pd(i,j) = sum(P(:,i)*Pd(:,j));
        end
    end
end

function dd = elemPdPd(Pd)
    % elemPdPd - integral ∫P'*P'dy of basis functions P (element flux)
    N = size(Pd,2);
    dd = zeros(N);
    for i = 1:N
        for j = 1:N
            dd(i,j) = -sum(Pd(:,i)*Pd(:,j));
            dd(j,i) = dd(i,j);
        end
    end
end

