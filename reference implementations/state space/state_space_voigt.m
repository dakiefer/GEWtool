%% Dispersion calculation in a generally anisotropic elastic plate
%% Example
% Compute dispersion spectrum for a 1-mm thick brass plate:

mat = Material('brass'); c = mat.tensor; rho = mat.rho;
h = 1e-3; % thickness in m
freq = linspace(0.001, 4, 500).'*1e6; % frequency in Hz
% freq = linspace(0.005, 4, 150)*1e6; % frequency in Hz
N = 15; % discretization = polynomial order of interpolants
%% 
% Define normalization parameters and convert physical quantities to normalized 
% quantities:

tic
c0 = c(1,2,1,2); h0 = h; % f0 = 1/2*(freq(1) + freq(end)); rho0 = c0/(f0^2*h0^2); % normalization parameters
rho0 = rho; f0 = sqrt(c0/rho0)/h0; % normalization parameters
wn = 2*pi*freq/f0; 
hn = h/h0;
rhon = rho/rho0; 
cn = c/c0;
polLamb = [1, 2, 6];
Cn = mat.C(polLamb, polLamb)/c0;

%% 
% Wave operator matrices $W = (i k_x h)^2 L_2 + i k_x h \partial_y L_1 + \partial_y^2 
% L_0 + (w_n h_n)^2 M$ and traction operator $\mathcal{T} = i k_x B_1 + \partial_y 
% B_0$

Iv = eye(2); IT = eye(3);
B1 = [1, 0; 0, 0; 0, 1];
B0 = [0, 0; 0, 1; 1, 0];
L0 = [0*Iv, B0.'; B0, 0*IT];
L1 = [0*Iv, B1.'; B1, 0*IT];
% L1 = [0*Iv, B1.'; 
%       Cn*B1, IT];
% L0 = [0*Iv, B0.';
%       Cn*B0, IT];
M = blkdiag(rhon*Iv, inv(Cn));
% BC = diag([0, 0, 0, 0, 1]);

% L0 = inv(M)*L0;
% L1 = inv(M)*L1;
% M = eye(size(M));

%% 
% Compute differentiation matrices on normalized domain $[h_0/2, h_0/2]$

[y_dash, Dy_dash] = chebdif(N, 2);
y = -h/2*y_dash;
Dy = -2*Dy_dash(:,:,1); % differentiation on [-h/2, h/2]
% e = ones(N,1);
% Dy = full(spdiags([-1/2*e, 0*e, 1/2*e], -1:1, N, N));
Id = eye(size(Dy));  % identity matrix for discretization
%% 
% Discretize and incorporate boundary conditions (traction free boundaries):
L1d = kron(L1,Id); L0d = kron(L0,-1i*Dy); Md = kron(M,Id);
% dofBC = [1, N];
% Bd = kron(B,Id(dofBC, :)); 

% dofBC = 3*N + [1, N, N+1, 2*N, 2*N+1, 3*N];
dofBC = 3*N + [1, N, N+1, 2*N];
L1d(dofBC, :) = []; L0d(dofBC, :) = []; Md(dofBC, :) = []; % remove DBC rows
L1d(:, dofBC) = []; L0d(:, dofBC) = []; Md(:, dofBC) = []; % remove DBC columns
% L1d(dofBC, :) = 0; L0d(dofBC, :) = Bd; Md(dofBC, :) = 0;
%% 
% Iterate over frequencies, solve and plot phase velocity

kh = nan(length(freq), size(Md, 2));
for ii = 1:length(freq)
    whn = wn(ii)*hn; % current frequency-thickness (normalized)
    [un, khi] = polyeig(L0d + whn*Md, L1d); % solve eigenvalue problem
    kh(ii, :) = khi; % extract kh
end
chron = toc;
fprintf('nF: %d, nK: %d, elapsed time: %g, time per point: %g.\n', size(kh, 1), size(kh, 2), chron, chron/(size(kh, 1)*size(kh, 2)));

%% 
% Prepare for plotting only propagating modes:

ff = repmat(freq, 1, size(kh, 2));
ind = abs(imag(kh)) < 1e-2; % indices of wavenumbers with low imaginary part
khp = kh(ind); fp = ff(ind);
%% 
% Plot real part of wavenumbers:

% figure, hold on, title('wavenumbers (real)'), xlim([0, 12]) 
% plot(real(khp), fp*h, 'k.')
% xlabel('Re{kx}*h in rad'), ylabel('f h in MHz mm')
%% 
% Plot phase velocity:

hold on, title('phase vel'), ylim([0, 10])
plot(fp*h, 2*pi*fp*h./real(khp)/1e3, 'r.')
yline(mat.cl/1e3), yline(mat.ct/1e3)
xlabel('f h in MHz mm'), ylabel('cp in mm/us')
