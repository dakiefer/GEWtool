%% compute axially guided ultrasonic waves in cylindrical coordinates
% 
% 2022 Daniel Kiefer
% Institut Langevin, Paris, France
% 
% see also:
% [1] M. Zheng, C. He, Y. Lyu, and B. Wu, “Guided waves propagation in anisotropic 
% hollow cylinders by Legendre polynomial solution based on state-vector formalism,” 
% Composite Structures, vol. 207, pp. 645–657, Jan. 2019, 
% doi: 10.1016/j.compstruct.2018.09.042.
% 

% specify parameters:
mat = jsondecode(fileread('../../../Material/database/steel_zheng.json'));
rho = mat.rho; c = voigt2tensor(mat.C);
% mat = Material('steel');
% rho = mat.rho; c = mat.tensor;
b = 141.3/2*1e-3; % outer radius
h = 12.5e-3; % thickness
a = b-h; % inner radius
n = 0; % circumferential wavenumber (flexural order)
N = 20; % number of collocation points
% fh = linspace(1e-2, 6, 300).'*1e3;
kh = linspace(1e-2, 15, 200);

%% derived parameters and normalize parameters:
h = b-a; % thickness
c0 = c(1,2,1,2); h0 = h; % normalization parameters
rho0 = rho; f0 = sqrt(c0/rho0)/h0; % normalization parameters
rhon = rho/rho0; cn = c/c0;

% relevant material matrices: 
crr = squeeze(cn(1,:,:,1));
cpp = squeeze(cn(2,:,:,2));
czz = squeeze(cn(3,:,:,3));
Czp = squeeze(cn(3,:,:,2)) + squeeze(cn(2,:,:,3));
Crz = squeeze(cn(3,:,:,1)) + squeeze(cn(1,:,:,3));
Crp = squeeze(cn(1,:,:,2)) + squeeze(cn(2,:,:,1));
crp = squeeze(cn(1,:,:,2));
crz = squeeze(cn(1,:,:,3));
A = [0, -1, 0; 1, 0, 0; 0, 0, 0]; % differetiation in curvilinear coordinate system
B = [1,  0, 0; 0, 1, 0; 0, 0, 0]; % differetiation in curvilinear coordinate system
rho = rhon*eye(size(crr));

%% discretization 
[y_dash, Dr_dash] = chebdif(N, 2); % y_dash in [-1, 1]
r = (h*y_dash + b + a)/2;
rn = r/h0;
Dr1 =  2*h0/h*Dr_dash(:,:,1); % differentiation on [-1/2, 1/2]
Dr2 = (2*h0/h)^2*Dr_dash(:,:,2); % 2nd order differentiation on [-1/2, 1/2]
Id = eye(size(Dr1));  % identity matrix for discretization
rn1inv = diag(1./rn);

%% problem setup
L2 = kron(czz, Id); 
L1 = kron(Czp*A, diag(1./rn)) + kron(Crz, Dr1) + 1i*n*kron(Czp, diag(1./rn));
L0 = kron(crr, Dr2) + kron(Crp*A, diag(1./rn)*Dr1 - diag(1./rn.^2)) ...
    + kron(cpp, diag(1./rn)*Dr1) - kron(cpp*B, diag(1./rn.^2)) ...
    + 1i*n*(kron(Crp, diag(1./rn)*Dr1 - diag(1./rn.^2)) + kron(2*cpp*A, diag(1./rn.^2))) ...
    + (1i*n).^2*(kron(cpp, diag(1./rn.^2)));
M = kron(rho, Id);
B1 = kron(crz, Id([1, N], :)); % BC
B0 = kron(crr, Dr1([1, N], :)) + kron(crp*A, rn1inv([1, N], :)) + 1i*n*kron(crp, rn1inv([1, N], :)); % BC going into L0

% incorporate BCs
dofBC = [1, N, N+1, 2*N, 2*N+1, 3*N];
L2(dofBC, :) = 0; L1(dofBC, :) = B1; L0(dofBC, :) = B0; M(dofBC, :) = 0;


%% solve for frequency:
tic 
nK = size(M, 2);
whn = nan(nK, length(kh));
for ii = 1:length(kh)
    kh0 = kh(ii);
    [wh2] = polyeig((1i*kh0)^2*L2 + (1i*kh0)*L1 + L0, M); % does not work properly with eig()
    whn(:,ii) = sqrt(wh2);
end
fh = real(whn/2/pi*f0*h0); fh(fh == 0) = nan;
chron = toc;
fprintf('nF: %d, nK: %d, elapsed time: %g, time per point: %g. ms\n', size(fh, 2), size(fh, 1), chron, chron/length(fh(:))*1e3);

% wavenumbers:
% kkh = kh.*ones(size(fh));
% figure, plot(kkh(:), fh(:), '.');
% ylim([0, 6e3]),
% ylabel('fh in m/s'), xlabel('kh in 1')

load('zheng.mat');
figure, hold on
plot(L01.f, L01.cp, 'k'), plot(L02.f, L02.cp, 'k'), plot(L04.f, L04.cp, 'k')
plot(F11.f, F11.cp, 'k'), plot(F12.f, F12.cp, 'k'), plot(F13.f, F13.cp, 'k')

% phase vel:
kkh = kh.*ones(size(fh));
plot(fh(:)/h0/1e6, 2*pi*fh(:)./kkh(:), 'r.');
ylim([0, 12e3]), xlim([0, 0.4])
xlabel('f in MHz'), ylabel('cp in m/s')

%% flexural waves of order 1
n = 1; % circumferential wavenumber (flexural order)


%% problem setup: (i*kh)^2*L2 + (i*kh)*L1(i*n) + L0(i*n) + w^2*M = 0
L2 = kron(czz, Id); 
L1 = kron(Czp*A, diag(1./rn)) + kron(Crz, Dr1) + 1i*n*kron(Czp, diag(1./rn));
L0 = kron(crr, Dr2) + kron(Crp*A, diag(1./rn)*Dr1 - diag(1./rn.^2)) ...
    + kron(cpp, diag(1./rn)*Dr1) - kron(cpp*B, diag(1./rn.^2)) ...
    + 1i*n*(kron(Crp, diag(1./rn)*Dr1 - diag(1./rn.^2)) + kron(2*cpp*A, diag(1./rn.^2))) ...
    + (1i*n).^2*(kron(cpp, diag(1./rn.^2)));
M = kron(rho, Id);

% incorporate BCs:
B1 = kron(crz, Id([1, N], :)); % BC going into L1
B0 = kron(crr, Dr1([1, N], :)) + kron(crp*A, rn1inv([1, N], :)) + 1i*n*kron(crp, rn1inv([1, N], :)); % BC going into L0
dofBC = [1, N, N+1, 2*N, 2*N+1, 3*N];
L2(dofBC, :) = 0; L1(dofBC, :) = B1; L0(dofBC, :) = B0; M(dofBC, :) = 0;

%% solve
tic 
nK = size(M, 2);
whn = nan(nK, length(kh));
for ii = 1:length(kh)
    kh0 = kh(ii);
    [wh2] = polyeig((1i*kh0)^2*L2 + (1i*kh0)*L1 + L0, M); % does not work properly with eig()
    whn(:,ii) = sqrt(wh2);
end
fh = real(whn/2/pi*f0*h0); fh(fh == 0) = nan;
chron = toc;
fprintf('nF: %d, nK: %d, elapsed time: %g, time per point: %g. ms\n', size(fh, 2), size(fh, 1), chron, chron/length(fh(:))*1e3);

% phase vel:
kkh = kh.*ones(size(fh));
plot(fh(:)/h0/1e6, 2*pi*fh(:)./kkh(:), 'g.');
ylim([0, 12e3]), xlim([0, 0.4])
xlabel('f in MHz'), ylabel('cp in m/s')
legend({'L01', 'L02', 'L04', 'F11', 'F12', 'F13', 'n=0', 'n=1'})
