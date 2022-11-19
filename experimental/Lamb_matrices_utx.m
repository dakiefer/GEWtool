function [gew] = Lamb_matrices_utx(mat, h, N)
% Generate matrices for convergence test
c = mat.c; rho = mat.rho;
c0 = c(1,2,1,2); h0 = h; % normalization parameters
rho0 = rho; fh0 = sqrt(c0/rho0); % normalization parameters
rhon = rho/rho0; cn = c/c0;

% relevant material matrices: 
udof = [1, 2];
cxx = squeeze(cn(1,udof,udof,1));
cxy = squeeze(cn(1,udof,udof,2)); 
cyx = squeeze(cn(2,udof,udof,1));
cyy = squeeze(cn(2,udof,udof,2));
I = eye(size(cxx));

%% discretization 
[y_dash, Dy_dash] = chebdif(N, 2); y = -h/2*y_dash;
Dy1 = -2*Dy_dash(:,:,1); % differentiation on [-h/2, h/2]
Dy2 = 4*Dy_dash(:,:,2); % differentiation on [-h/2, h/2]
Id = eye(size(Dy1));  % identity matrix for discretization

%% problem setup: (i*kh)^2*L2 + (i*kh)*L1(i*n) + L0(i*n) + w^2*M = 0
%  (u, tx): displacements, axial tractions
L1 = [ kron(cyx, Dy1), kron(h0*I, Id); 
       kron(cxx, Id),  kron(0*I, Id)];
L0 = [ kron(cyy, Dy2), kron(0*I, Id); 
       kron(cxy, Dy1), kron(-h0*I, Id)];
M = blkdiag(kron(rhon*I, Id), kron(0*I, Id));

B1 = [kron(cyx, Id([1, N], :)), kron(0*I, Id([1, N], :))];
B0 = [kron(cyy, Dy1([1, N], :)), kron(0*I, Id([1, N], :))];
dofBC = [1, N, N+1, 2*N];
L1(dofBC, :) = B1; L0(dofBC, :) = B0; M(dofBC, :) = 0; 
% remove unused dofs (this is not necessary but reduces the matrix size):
dof = setdiff(1:4*N, [2*N+1, 3*N, 3*N+1, 4*N]);
L1 = L1(dof, dof); L0 = L0(dof, dof); M = M(dof, dof);

L2 = zeros(size(M)); 

gew.op.L0 = L0; 
gew.op.L1 = L1;
gew.op.L2 = L2;
gew.op.M = M;
gew.np.fh0 = fh0;
gew.np.h0 = h0;
gew.np.c0 = c0;
gew.np.rho0 = rho0;
gew.geom = Geometry([0, h],N,2);

end
