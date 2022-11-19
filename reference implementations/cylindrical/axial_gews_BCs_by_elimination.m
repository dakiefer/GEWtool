%% compute axially guided ultrasonic waves in cylindrical coordinates
% Impose BCs by elimination, see [2]. This leads to regular matrices L and M. 
% For some reason, eig() is still not always able to find the solutions.
% 
% 2022 Daniel Kiefer
% Institut Langevin, Paris, France
% 
% [2] J. Hoepﬀner, “Implementation of boundary conditions,” Paris, 2007. [Online]. 
% Available: http://www.lmm.jussieu.fr/~hoepffner/home.php
% 

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
rho = rhon*eye(size(crr)); % expand to matrix

%% discretization 
[y_dash, Dr_dash] = chebdif(N, 2);  % y_dash in [-1, 1]
r = (h*y_dash + b + a)/2;           % radial coordinate
rn = r/h0;
Dr1 =  2*h0/h*Dr_dash(:,:,1);       % differentiation on [-1/2, 1/2]
Dr2 = (2*h0/h)^2*Dr_dash(:,:,2);    % 2nd order differentiation on [-1/2, 1/2]
Id = eye(size(Dr1));                % identity matrix for discretization
rn1inv = diag(1./rn);               % 1/r^2 (mostly for use in BC)

%% problem setup: (i*kh)^2*L2 + (i*kh)*L1(i*n) + L0(i*n) + w^2*M = 0
L2 = kron(czz, Id); 
L1 = kron(Czp*A, diag(1./rn)) + kron(Crz, Dr1) + 1i*n*kron(Czp, diag(1./rn));
L0 = kron(crr, Dr2) + kron(Crp*A, diag(1./rn)*Dr1 - diag(1./rn.^2)) ...
    + kron(cpp, diag(1./rn)*Dr1) - kron(cpp*B, diag(1./rn.^2)) ...
    + 1i*n*(kron(Crp, diag(1./rn)*Dr1 - diag(1./rn.^2)) + kron(2*cpp*A, diag(1./rn.^2))) ...
    + (1i*n).^2*(kron(cpp, diag(1./rn.^2)));
M = kron(rho, Id);

% define BCs:
B1 = kron(crz, Id([1, N], :)); % BC going into L1
B0 = kron(crr, Dr1([1, N], :)) + kron(crp*A, rn1inv([1, N], :)); % BC going into L0
rem = [1, N, N+1, 2*N, 2*N+1, 3*N]; % dofs to remove
dof = setdiff(1:size(M,1), rem);    % dofs to keep
G = @(ikh) -(ikh*B1(:,rem) + B0(:,rem))\(ikh*B1(:,dof) + B0(:,dof)); % give-back matrix: u_dofBC = G*u_dof

%% solve for frequency:
tic 
whn = nan(length(dof), length(kh));
for ii = 1:length(kh)
    khi = kh(ii);
    Li = (1i*khi)^2*L2 + (1i*khi)*L1 + L0;
    Li = Li(dof, dof) + Li(dof, rem)*G(1i*khi); % incorporate BC
    Mi = M(dof, dof) + M(dof, rem)*G(1i*khi);
    [wh2] = eig(Li, Mi); % does not work properly with eig()
    whn(:,ii) = -1i*sqrt(wh2);
    if cond(Mi) > 1e6 || cond(Li) > 1e6, warning('close to singular'); end
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
% kkh = kh.*ones(size(fh));
% figure, plot(fh(:), 2*pi*fh(:)./kkh(:), '.');
% ylim([0, 15e3]), xlim([0, 6e3])
% xlabel('fh in m/s'), ylabel('cp in m/s')
