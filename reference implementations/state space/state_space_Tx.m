%% Dispersion calculation in a generally anisotropic elastic plate

% parameters:
mat = Material('brass'); c = mat.tensor; rho = mat.rho;
h = 1e-3; % thickness in m
kh = linspace(1e-2, 20, 200); % wavenumber-thickness (solve for frequency)
N = 20; % discretization = polynomial order of interpolants

%% derived parameters and normalize parameters:
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

%% displacements only
% L2 = kron(cxx, Id); L1 = kron((cxy + cyx), Dy1); L0 = kron(cyy, Dy2); M = kron(rhon*I, Id);
% B1 = kron(cyx, Id([1, N], :)); B0 = kron(cyy, Dy1([1, N], :)); 
% dofBC = [1, N, N+1, 2*N];
% L2(dofBC, :) = 0; L1(dofBC, :) = B1; L0(dofBC, :) = B0; M(dofBC, :) = 0;

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
L2 = zeros(size(L1)); 
% remove unused dofs (this is not necessary):
% dof = setdiff(1:4*N, [2*N+1, 3*N, 3*N+1, 4*N]);
% L1 = L1(dof, dof); L0 = L0(dof, dof); M = M(dof, dof); L2 = L2(dof, dof);

%% solve for frequency:
tic 
whn = nan(size(M, 2), length(kh));
for ii = 1:length(kh)
    kh0 = kh(ii);
    [wh2] = polyeig((1i*kh0)^2*L2 + (1i*kh0)*L1 + L0, M); % does not work properly with eig()
    whn(:,ii) = sqrt(wh2);
end
fh = real(whn/2/pi*fh0); % fh(fh == 0) = nan;
chron = toc;
fprintf('nF: %d, nK: %d, elapsed time: %g, time per point: %g. ms\n', size(fh, 2), size(fh, 1), chron, chron/length(fh(:))*1e3);

% wavenumbers:
% kkh = kh.*ones(size(fh));
% hold on, plot(kkh(:), fh(:), '.');
% ylim([0, 5e3]), xlim([0, 10])
% ylabel('fh in m/s'), xlabel('kh in 1')

% phase vel:
kkh = kh.*ones(size(fh));
figure, plot(fh(:), 2*pi*fh(:)./kkh(:), '.');
ylim([0, 10e3]), xlim([0, 6000])
xlabel('f h in m/s'), ylabel('cp in m/s')
