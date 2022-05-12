function [guw] = Lamb_matrices_rectangularSCM(mat, h, N)

%% Eigenvalue problem describing the ZGV points 
% Generates the matrices describing the Lamb wave problem as well as the associated 
% equation describing the ZGV points [1].
% Aluminum: S1S2: (kh = 1.62, fh = 2842 m/s), S3S6: (kh = 1.72, fh = 9311 m/s)
% 
% 2022 - Daniel Kiefer
% Institut Langevin, Paris, France
% 
% [1] J. L. Tassoulas and T. R. Akylas, “On Wave Modes With Zero Group Velocity in 
%     an Elastic Layer,” Journal of Applied Mechanics, vol. 51, no. 3, pp. 652–656, 
%     Sep. 1984, doi: 10.1115/1.3167688.

% parameters: 
c = mat.c; rho = mat.rho;

%% normalize parameters 
c0 = c(1,2,1,2); h0 = h; % normalization parameters
rho0 = rho; f0 = sqrt(c0/rho0)/h0; % normalization parameters
rhon = rho/rho0; cn = c/c0;

% relevant stiffness tensors: 
cxx = squeeze(cn(1,1:2,1:2,1));
cxy = squeeze(cn(1,1:2,1:2,2));
cyx = squeeze(cn(2,1:2,1:2,1));
cyy = squeeze(cn(2,1:2,1:2,2));
Rho = rhon*eye(size(cxx));

%% discretization 
% [~, Dy_dash] = chebdif(N, 2);
% Dy1 = -2*Dy_dash(:,:,1); % differentiation on [-1/2, 1/2]
% Dy2 = (2)^2*Dy_dash(:,:,2); % 2nd order differentiation on [-1/2, 1/2]
Dy1 = diffmat(N, 1, [-0.5 0.5], 'chebkind2'); % first order on domain [0 1]
Dy2 = diffmat(N, 2, [-0.5 0.5], 'chebkind2'); % second order on domain [0 1]
Id = eye(size(Dy1));  % identity matrix for discretization
[xf, ccw, baryw] = chebpts(N, [-0.5 0.5], 2); % second-kind points (includes boundaries)
[xd] = chebpts(N-2, [-0.5 0.5], 1); % first-kind points (without boundaries)
P = barymat(xd, xf, baryw); % resampling matrix (from values on xf to values on xd)
% P = eye(size(Id)); % disable resampling

%% Lamb wave problem
L2 = kron(cxx, P*Id); L1 = kron(cxy + cyx, P*Dy1); L0 = kron(cyy, P*Dy2); M = kron(Rho,P*Id);
B1 = kron(cyx, Id([1, N], :)); B0 = kron(cyy, Dy1([1, N], :)); % BCs
% dofBC = [1, N, N+1, 2*N];
% L2(dofBC, :) = B1; L1(dofBC, :) = B0; L0(dofBC, :) = 0; M(dofBC, :) = 0;
L2 = [zeros(size(B0)); L2]; L1 = [B1; L1]; L0 = [B0; L0]; M = [zeros(size(B0)); M]; % append BCs

guw.op.L0 = L0; 
guw.op.L1 = L1;
guw.op.L2 = L2;
guw.op.M = M;
guw.np.fh0 = f0*h0;
guw.np.h0 = h0;
guw.np.c0 = c0;
guw.np.rho0 = rho0;
guw.geom = Geometry([0, h],N,2);
