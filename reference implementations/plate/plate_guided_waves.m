%% Dispersion calculation in a generally anisotropic elastic plate
% Depends on the DMSUITE toolbox by Weideman and Reddy: 
% https://mathworks.com/matlabcentral/fileexchange/29-dmsuite
% 
% 2022 - Daniel A. Kiefer
% Institut Langevin, Paris, France
% 
% see also:
% [1] D. A. Kiefer, M. Ponschab, S. J. Rupitsch, and M. Mayle, “Calculating the full 
% leaky Lamb wave spectrum with exact fluid interaction,” The Journal of the Acoustical 
% Society of America, vol. 145, no. 6, pp. 3341–3350, Jun. 2019, doi: 10.1121/1.5109399.
%

h = 1e-3;   % thickness in m
N = 20;     % discretization: polynomial order of interpolants
rho = 7900; lbd = 1.1538e11; mu = 7.6923e10; % steel material

% define and normalize parameters:
II = eye(3).*shiftdim(eye(3), -2); % 4th order "unit tensor"
c = lbd*II + mu*(permute(II, [1 3 4 2]) + permute(II, [1 3 2 4])); % stiffness tensor
c0 = c(1,2,1,2); h0 = h; % normalization parameters
rho0 = rho; fh0 = sqrt(c0/rho0); % normalization parameters
rhon = rho/rho0; cn = c/c0;

% relevant material matrices: 
udof = 1:3; % Lamb and/or SH
cxx = squeeze(cn(1,udof,udof,1));
cxy = squeeze(cn(1,udof,udof,2)); 
cyx = squeeze(cn(2,udof,udof,1));
cyy = squeeze(cn(2,udof,udof,2));
I = eye(size(cxx)); 

%% discretization 
tic
[~, Dy_dash] = chebdif(N, 2); 
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
chron = toc; fprintf('assembling time: %g\n', chron);

%% solve for frequency and plot:
kh = linspace(1e-2, 15, 300); % wavenumber*thickness 
whn = nan(size(M, 2), length(kh)); tic 
for ii = 1:length(kh)
    kh0 = kh(ii);
    [wh2] = polyeig((1i*kh0)^2*L2 + (1i*kh0)*L1 + L0, M); 
    whn(:,ii) = sqrt(wh2);
end
fh = real(whn/2/pi*fh0); fh(fh == 0) = nan;
chron = toc; fprintf('nF: %d, nK: %d, elapsed time: %g, time per point: %g. ms\n', size(fh, 1), size(fh, 2), chron, chron/length(fh(:))*1e3);

% plot wavenumbers:
kkh = kh.*ones(size(fh));
figure(1), plot(kkh(:)/h/1e3, fh(:)/h/1e6, '.');
xlim([0, 12]), ylim([0, 6]),
xlabel('k in rad/mm'), ylabel('f in MHz'),

% plot phase vel:
kkh = kh.*ones(size(fh));
figure(2), plot(fh(:)/h/1e6, 2*pi*fh(:)./kkh(:)/1e3, '.');
xlim([0, 6e6]/1e6), ylim([0, 12])
xlabel('f in MHz'), ylabel('cp in mm/us')

%% solve for wavenumbers:
fh = linspace(1e-2, 6, 300).'*1e6*h; % frequency*thickness
kh = nan(length(fh), size(M, 2)*2); tic
for ii = 1:length(fh)
    whn = 2*pi*fh(ii)/fh0; % current frequency-thickness (normalized)
    [un, khi] = polyeig(L0 + whn^2*M, L1, L2); 
    kh(ii, :) = -1i*khi;
end
chron = toc; fprintf('nF: %d, nK: %d, elapsed time: %g, time per point: %g. ms\n', size(kh, 1), size(kh, 2), chron, chron/length(kh(:))*1e3);

% plot wave numbers:
ffh = fh.*ones(size(kh));
figure(3), scatter(real(kh(:))/h/1e3, ffh(:)/h/1e6, 8, abs(imag(kh(:))), 'filled'), 
caxis([0, 0.12]), xlim([0, 12]), ylim([0, fh(end)/h/1e6])
xlabel('k in rad/mm'), ylabel('f in MHz')

% plot phase vel:
ffh = fh.*ones(size(kh));
figure(4), scatter(ffh(:)/h/1e6, 2*pi*ffh(:)./real(kh(:))/1e3, 8, abs(imag(kh(:))), 'filled'), 
caxis([0, 0.12]), xlim([0, fh(end)/h/1e6]), ylim([0, 15]),
ylabel('cp in mm/us'), xlabel('f in MHz')

