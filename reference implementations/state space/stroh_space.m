%% Dispersion calculation in a generally anisotropic elastic plate
%% Example
% Compute dispersion spectrum for a 1-mm thick brass plate:

mat = Material('brass'); c = mat.c; rho = mat.rho;
h = 1e-3; % thickness in m
vp = linspace(0.001, 25, 500).'*1e3; % frequency in Hz
% freq = linspace(0.005, 4, 150)*1e6; % frequency in Hz
N = 15; % discretization = polynomial order of interpolants
%% 
% Define normalization parameters and convert physical quantities to normalized 
% quantities:

tic
c0 = c(1,2,1,2); h0 = h; % f0 = 1/2*(freq(1) + freq(end)); rho0 = c0/(f0^2*h0^2); % normalization parameters
rho0 = rho; f0 = sqrt(c0/rho0)/h0; % normalization parameters
vp0 = h0*f0;
vpn = vp/vp0;
hn = h/h0;
rhon = rho/rho0; 
cn = c/c0;
%% 
% Wave operator matrices $W = (i k_x h)^2 L_2 + i k_x h \partial_y L_1 + \partial_y^2 
% L_0 + (w_n h_n)^2 M$ and traction operator $\mathcal{T} = i k_x B_1 + \partial_y 
% B_0$

xy = [1, 2];
cxx = squeeze(cn(1,xy,xy,1));
cxy = squeeze(cn(1,xy,xy,2)); cyx = squeeze(cn(2,xy,xy,1));
cyy = squeeze(cn(2,xy,xy,2));
I = eye(size(cxx)); 
Bvp = [rhon*I,  0*I;
       0*I,     0*I];
B0 = [-cxx, 0*I;
      -cyx,  I];
A = [cxy,   I  ; 
     cyy,   0*I];

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
Ad = kron(A,Dy); B0d = kron(B0,Id); Bvpd = kron(Bvp,Id);
dofBC = 2*N + [1, N, N+1, 2*N];
Ad(dofBC, :) = []; B0d(dofBC, :) = []; Bvpd(dofBC, :) = []; % remove DBC rows
Ad(:, dofBC) = []; B0d(:, dofBC) = []; Bvpd(:, dofBC) = []; % remove DBC columns
% L1d(dofBC, :) = 0; L0d(dofBC, :) = Bd; Md(dofBC, :) = 0;
%% 
% Iterate over frequencies, solve and plot phase velocity

kh = nan(length(vp), size(Ad, 2));
for ii = 1:length(vpn)
    vpi = vpn(ii); % current phase velocity
    [V, D] = eig(Ad, vpi^2*Bvpd+B0d); % solve eigenvalue problem
    khi = -1i*diag(D);
    kh(ii, :) = khi; % extract kh
end
chron = toc;
fprintf('nF: %d, nK: %d, elapsed time: %g, time per point: %g.\n', size(kh, 1), size(kh, 2), chron, chron/(size(kh, 1)*size(kh, 2)));
%% 
% Prepare for plotting only propagating modes:

ff = real(vpn*vp0.*kh/h)/2/pi;
ind = abs(imag(kh)) < 1e-2; % indices of wavenumbers with low imaginary part
khp = kh(ind); fp = ff(ind);
vpp = repmat(vp, 1, size(kh, 2));
vpp = vpp(ind);
%% 
% Plot real part of wavenumbers:

figure, hold on, title('wavenumbers (real)')
plot(real(khp/h), fp, 'k.')
xlabel('Re{kx} in rad/m'), ylabel('f in MHz')
xlim([0, 12e3]), ylim([0, 6e6])
%% 
% Plot phase velocity:
figure, hold on, title('phase vel'), %ylim([0, 10])
plot(fp, vpp/1e3, 'k.')
yline(mat.cl/1e3), yline(mat.ct/1e3)
xlabel('f in MHz'), ylabel('cp in mm/us')
xlim([0, 6e6])
