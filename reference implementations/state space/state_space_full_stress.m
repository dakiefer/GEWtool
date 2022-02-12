%% Dispersion calculation in a generally anisotropic elastic plate
%% Example
% Compute dispersion spectrum for a 1-mm thick brass plate:

mat = Material('brass'); c = mat.tensor; rho = mat.rho;
h = 1e-3; % thickness in m
freq = linspace(0.001, 4, 500).'*1e6; % frequency in Hz
% freq = 1e6;
% freq = linspace(0.005, 4, 150)*1e6; % frequency in Hz
N = 20; % discretization = polynomial order of interpolants
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
%% 
% Wave operator matrices $W = (i k_x h)^2 L_2 + i k_x h \partial_y L_1 + \partial_y^2 
% L_0 + (w_n h_n)^2 M$ and traction operator $\mathcal{T} = i k_x B_1 + \partial_y 
% B_0$

xy = [1, 2];
cxx = squeeze(cn(1,xy,xy,1));
cxy = squeeze(cn(1,xy,xy,2)); cyx = squeeze(cn(2,xy,xy,1));
cyy = squeeze(cn(2,xy,xy,2));
I = eye(size(cxx)); I2 = [0, 0; 1, 0];
L1 = [0*I,  1/rhon*I, I2; 
      cxx, 0*I, 0*I; 
      cyx, 0*I, 0*I];
L0 = -1i*[0*I,   0*I,   1/rhon*I; 
          cxy, 0*I,   0*I; 
          cyy, 0*I,   0*I];
M = [I,   0*I,   0*I; 
     0*I, I,     I2; 
     0*I, 0*I,   I];
B = [0*I, 0*I, I]; % boundary condition term
dofTxy = 4;
L1(dofTxy, :) = []; L0(dofTxy, :) = []; M(dofTxy, :) = [];
L1(:, dofTxy) = []; L0(:, dofTxy) = []; M(:, dofTxy) = [];
B(:, dofTxy) = [];
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
L1d = kron(L1,Id); L0d = kron(L0,Dy); Md = kron(M,Id);
dofBC = [1, N];
Bd = kron(B,Id(dofBC, :)); 

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
    [un, khiMat] = eig(L0d + whn*Md, L1d);
%     unorm = sum(un([1:2*N], :).*conj(un([1:2*N], :)), 1);
%     Tnorm = sum(un([2*N+1:size(un, 1)], :).*conj(un([2*N+1:size(un, 1)], :)), 1);
    khi = diag(khiMat);
%     khi(unorm./Tnorm < 1e-1) = nan;
%     [un, khi] = polyeig(L0d + whn*Md, L1d); % solve eigenvalue problem
    kh(ii, :) = khi; % extract kh
%     [~, ki] = min(2*pi*freq*h0./kh - mat.cl);
%     plot(abs(un(:, ki))); xline(2*N)
end
chron = toc;
fprintf('nF: %d, nK: %d, elapsed time: %g, time per point: %g.\n', size(kh, 1), size(kh, 2), chron, chron/(size(kh, 1)*size(kh, 2)));

% 
% Prepare for plotting only propagating modes:

ff = repmat(freq, 1, size(kh, 2));
ind = abs(imag(kh)) < 1e-2; % indices of wavenumbers with low imaginary part
khp = kh(ind); fp = ff(ind);

clf, hold on, title('phase vel'), ylim([0, 10])
plot(fp*h, 2*pi*fp*h./real(khp)/1e3, 'k.')
yline(mat.cl/1e3), yline(mat.ct/1e3)
xlabel('f h in MHz mm'), ylabel('cp in mm/us')
