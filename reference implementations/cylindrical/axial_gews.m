%% compute axially guided ultrasonic waves in cylindrical coordinates
% 
% see also:
% [1] M. Zheng, C. He, Y. Lyu, and B. Wu, “Guided waves propagation in anisotropic 
% hollow cylinders by Legendre polynomial solution based on state-vector formalism,” 
% Composite Structures, vol. 207, pp. 645–657, Jan. 2019, 
% doi: 10.1016/j.compstruct.2018.09.042.
% 
% 2022 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France

% specify parameters:
a = 3e-3; b = 5e-3; % inner and outer radius
n = 0; % circumferential wavenumber (flexural order)
N = 20; % number of collocation points
kh = linspace(1e-2, 15, 400); % wavenumber-thickness (solve for frequency)
% fh = linspace(1e-2, 6, 300).'*1e3; % frequency-thickness (solve for wavenumbers)
rho = 7900; lbd = 1.1538e11; mu = 7.6923e10; % material parameters
II = eye(3).*shiftdim(eye(3), -2); % 4th order "unit tensor"
c = lbd*II + mu*(permute(II, [1 3 4 2]) + permute(II, [1 3 2 4])); % stiffness tensor

%% derived parameters and normalize parameters:
h = b-a; % thickness
c0 = c(1,2,1,2); h0 = h; % normalization parameters
rho0 = rho; f0 = sqrt(c0/rho0)/h0; % normalization parameters
rhon = rho/rho0; cn = c/c0;
udof = [1, 2, 3];

% relevant material matrices: 
cxx = squeeze(cn(1,udof,udof,1));
crr = squeeze(cn(2,udof,udof,2));
cpp = squeeze(cn(3,udof,udof,3));
Crp = squeeze(cn(2,udof,udof,3)) + squeeze(cn(3,udof,udof,2));
Cxp = squeeze(cn(1,udof,udof,3)) + squeeze(cn(3,udof,udof,1));
Crx = squeeze(cn(2,udof,udof,1)) + squeeze(cn(1,udof,udof,2));
crp = squeeze(cn(2,udof,udof,3));
crx = squeeze(cn(2,udof,udof,1));
A = [0, 0, 0; 0, 0, -1; 0, 1, 0]; % differetiation in curvilinear coordinate system
B = [0, 0, 0; 0, 1,  0; 0, 0, 1]; % differetiation in curvilinear coordinate system
A = A(udof, udof); B = B(udof, udof);
rho = rhon*eye(size(crr)); % expand to matrix

% cpp = 0*cpp; Czp = 0*Czp; Crp = 0*Crp; crp = 0*crp;

%% discretization 
[y_dash, Dr_dash] = chebdif(N, 2);  % y_dash in [-1, 1]
r = (h*y_dash + b + a)/2;           % radial coordinate
rn = r/h;
Dr1 =  2*h0/h*Dr_dash(:,:,1);       % differentiation on [-1/2, 1/2]
Dr2 = (2*h0/h)^2*Dr_dash(:,:,2);    % 2nd order differentiation on [-1/2, 1/2]
Id = eye(size(Dr1));                % identity matrix for discretization
rn1inv = diag(1./rn);               % 1/r^2 (mostly for use in BC)

%% problem setup: (i*kh)^2*L2 + (i*kh)*L1(i*n) + L0(i*n) + w^2*M = 0
L2 = kron(cxx, Id); 
L1 = kron(Cxp*A, diag(1./rn)) + kron(Crx, Dr1) + 1i*n*kron(Cxp, diag(1./rn));
L0 = kron(crr, Dr2) + kron(Crp*A, diag(1./rn)*Dr1 - diag(1./rn.^2)) ...
    + kron(cpp, diag(1./rn)*Dr1) - kron(cpp*B, diag(1./rn.^2)) ...
    + 1i*n*(kron(Crp, diag(1./rn)*Dr1 - diag(1./rn.^2)) + kron(2*cpp*A, diag(1./rn.^2))) ...
    + (1i*n).^2*(kron(cpp, diag(1./rn.^2)));
M = kron(rho, Id);

% incorporate BCs:
B1 = kron(crx, Id([1, N], :)); % BC going into L1
B0 = kron(crr, Dr1([1, N], :)) + kron(crp*A, rn1inv([1, N], :)) + 1i*n*kron(crp, rn1inv([1, N], :)); % BC going into L0
dofBC = [(0:length(udof)-1)*N+1; (1:length(udof))*N]; % [1, N, N+1, 2*N, 2*N+1, 3*N];
L2(dofBC, :) = 0; L1(dofBC, :) = B1; L0(dofBC, :) = B0; M(dofBC, :) = 0;

%% solve for frequency:
tic 
whn = nan(size(M, 2), length(kh));
for ii = 1:length(kh)
    kh0 = kh(ii);
    [wh2] = polyeig((1i*kh0)^2*L2 + (1i*kh0)*L1 + L0, M); % does not work properly with eig()
    whn(:,ii) = sqrt(wh2);
end
fh = real(whn/2/pi*f0*h0); fh(fh == 0) = nan;
chron = toc;
fprintf('nF: %d, nK: %d, elapsed time: %g, time per point: %g. ms\n', size(fh, 2), size(fh, 1), chron, chron/length(fh(:))*1e3);

% wavenumbers:
kkh = kh.*ones(size(fh));
figure, plot(kkh(:), fh(:), '.');
ylim([0, 6e3]),
ylabel('fh in m/s'), xlabel('kh in 1')

% phase vel:
kkh = kh.*ones(size(fh));
figure, plot(fh(:), 2*pi*fh(:)./kkh(:), '.');
ylim([0, 15e3]), xlim([0, 6e3])
xlabel('fh in m/s'), ylabel('cp in m/s')

%% solve for wavenumbers:
% tic
% kh = nan(length(fh), size(M, 2)*2);
% for ii = 1:length(fh)
%     whn = 2*pi*fh(ii)/f0/h0; % current frequency-thickness (normalized)
%     [un, khi] = polyeig(L0 + whn^2*M, L1, L2); % solve eigenvalue problem
%     kh(ii, :) = -1i*khi; % extract kh
% end
% chron = toc;
% fprintf('nF: %d, nK: %d, elapsed time: %g, time per point: %g. ms\n', size(kh, 1), size(kh, 2), chron, chron/length(kh(:))*1e3);
% 
% % wave numbers:
% ffh = fh.*ones(size(kh));
% figure, scatter(real(kh(:)), ffh(:), 8, abs(imag(kh(:))), 'filled'), 
% caxis([0, 0.12]), ylim([0, fh(end)]), xlim([0, 14]),
% xlabel('kh in 1'), ylabel('fh in m/s')
% 
% % phase vel:
% ffh = fh.*ones(size(kh));
% figure, scatter(ffh(:), 2*pi*ffh(:)./real(kh(:)), 8, abs(imag(kh(:))), 'filled'), 
% caxis([0, 0.12]), xlim([0, fh(end)]), ylim([0, 15e3]),
% ylabel('cp in m/s'), xlabel('fh in m/s')

%% compare with reference implementation:
% plate = Plate(mat, h, 30);
% Nmodes = 80;
% freq = linspace(1e-2, 6, 300).'*1e6;
% data = plate.solveForFreq(freq, Nmodes);
% dr = data.selectModes(real(data.k) >= 0); 
% ff0 = freq.*ones(size(dr.k));
% plot(real(dr.k(:))*h, ff0(:)*h, 'rx'), hold on
% % scatter(real(dr.k(:)*h), ff(:)*h, 10, abs(imag(dr.k(:)*h)), 'filled'), caxis([0, 0.12])
% plot(kkh(:), fh(:), 'k.');
% ylim([0, 6e3]), xlim([0, kh(end)])
% ylabel('fh in m/s'), xlabel('kh in 1')
% legend({'plate', 'cylinder'}, 'Location', 'south east')