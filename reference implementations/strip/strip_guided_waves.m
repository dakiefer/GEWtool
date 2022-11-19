%% Dispersion calculation in a generally anisotropic elastic strip
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

h = 5e-3;   % thickness in m (ey-dimension)
b = 15e-3;   % width in m (ez-dimension)
N = 8;     % discretization in y: polynomial order of interpolants
P = 10;     % discretization in z: polynomial order of interpolants
rho = 7900; lbd = 1.1538e11; mu = 7.6923e10; % steel material

% define and normalize parameters:
II = eye(3).*shiftdim(eye(3), -2); % 4th order "unit tensor"
c = lbd*II + mu*(permute(II, [1 3 4 2]) + permute(II, [1 3 2 4])); % stiffness tensor
c0 = c(1,2,1,2); h0 = h; % normalization parameters
rho0 = rho; fh0 = sqrt(c0/rho0); % normalization parameters
rhon = rho/rho0; cn = c/c0;
hn = h/h0; bn = b/h0;

% relevant material matrices: 
udof = 1:3; % Lamb and/or SH
cxx = squeeze(cn(1,udof,udof,1));
cxy = squeeze(cn(1,udof,udof,2));
cxz = squeeze(cn(1,udof,udof,3));
cyx = squeeze(cn(2,udof,udof,1));
cyy = squeeze(cn(2,udof,udof,2));
cyz = squeeze(cn(2,udof,udof,3));
czx = squeeze(cn(3,udof,udof,1));
czy = squeeze(cn(3,udof,udof,2));
czz = squeeze(cn(3,udof,udof,3));
I = eye(size(cxx)); 

%% discretization 
[~, Dy_dash] = chebdif(N, 2); 
D10 = -2/hn*Dy_dash(:,:,1); % differentiation on normalized domain
D20 = 4/hn^2*Dy_dash(:,:,2);  % second order derivative
[~, Dz_dash] = chebdif(P, 2); 
D01 = -2/bn*Dz_dash(:,:,1); % differentiation on normalized domain
D02 = 4/bn^2*Dz_dash(:,:,2);  % second order derivative

% diff matrices on 2d-grid:
Dyz = kron(D01,D10); % = Dzy
Dy  = kron(eye(P),D10);
Dyy = kron(eye(P),D20);
Dz  = kron(D01,eye(N));
Dzz = kron(D02,eye(N));
Id =  eye(N*P);

% define wave operators:
L2 = kron(cxx, Id); 
L1 = kron(cyx + cxy, Dy) + kron(czx + cxz, Dz);
L0 = kron(cyy, Dyy) + kron(czy, Dyz) + kron(cyz, Dyz) + kron(czz, Dzz);
M  = kron(rhon*I, Id);
% define boundary operators:
% dofBC = kron([1 P], [1 N]);
Ind = reshape(1:length(Dy), N, P); % node enumeration
dofty = [Ind(1,:), Ind(N,:)]; % where to impose ty = ey.T = 0
doftz = [Ind(:,1); Ind(:,P)]; % where to impose tz = ez.T = 0
% boundary operators for traction ty:
By1 = kron(cyx, Id(dofty, :)); 
By0 = kron(cyy, Dy(dofty, :)) + kron(cyz, Dz(dofty, :));
% boundary operators for traction tz:
Bz1 = kron(czx, Id(doftz, :)); 
Bz0 = kron(czy, Dy(doftz, :)) + kron(czz, Dz(doftz, :));

% incorporate BCs:
doftyg = (N*P)*(0:length(udof)-1) + dofty.'; % dof for all three displ. ux,uy,uz
doftzg = (N*P)*(0:length(udof)-1) + doftz;
L2(doftyg, :) = 0; L1(doftyg, :) = By1; L0(doftyg, :) = By0; M(doftyg, :) = 0;
L2(doftzg, :) = 0; L1(doftzg, :) = Bz1; L0(doftzg, :) = Bz0; M(doftzg, :) = 0;

%% solve for frequency and plot:
kh = linspace(1e-2, 600, 120)*h0; % wavenumber*thickness 
whn = nan(size(M, 2), length(kh)); tic 
for ii = 1:length(kh)
    kh0 = kh(ii);
    [wh2] = polyeig((1i*kh0)^2*L2 + (1i*kh0)*L1 + L0, M); 
    whn(:,ii) = sqrt(wh2);
end
fh = real(whn/2/pi*fh0); fh(abs(fh) <= 1e-3) = nan;
chron = toc; fprintf('nF: %d, nK: %d, elapsed time: %g, time per point: %g. ms\n', size(fh, 1), size(fh, 2), chron, chron/length(fh(:))*1e3);

%% plot 
% for reference, plot plate first:
mat = Material('steel', c, rho);
plate = Plate(mat, h, 25);
gew = plate.fullyCoupled;
dat = computeW(gew, kh/h0);
figure, hold on
plot(dat.k(:), dat.w(:)/2/pi, 'kx');
plate = Plate(mat, b, 25);
gew = plate.fullyCoupled;
dat = computeW(gew, kh/h0);
plot(dat.k(:), dat.w(:)/2/pi, 'x', 'Color', [0.7 0.7 0.7]);

% plot strip wavenumbers:
kkh = kh.*ones(size(fh));
plot(kkh(:)/h0, fh(:)/h0, 'r.', 'MarkerSize', 10);
xlim([0, max(kh)/h0]), ylim([0, 300e3]),
xlabel('k in rad/m'), ylabel('f in Hz'),
legend({'plate h', 'plate b', 'strip hxb'}, 'Location','southeast')
title(sprintf('strip with cross section h = %d mm, b = %d mm', h/1e-3, b/1e-3))

% plot phase vel:
% kkh = kh.*ones(size(fh));
% figure(2), plot(fh(:)/h0/1e6, 2*pi*fh(:)./kkh(:)/1e3, '.');
% xlim([0, 6e6]/1e6), ylim([0, 12])
% xlabel('f in MHz'), ylabel('cp in mm/us')

% %% solve for wavenumbers:
% fh = linspace(1e-2, 6, 300).'*1e6*h; % frequency*thickness
% kh = nan(length(fh), size(M, 2)*2); tic
% for ii = 1:length(fh)
%     whn = 2*pi*fh(ii)/fh0; % current frequency-thickness (normalized)
%     [un, khi] = polyeig(L0 + whn^2*M, L1, L2); 
%     kh(ii, :) = -1i*khi;
% end
% chron = toc; fprintf('nF: %d, nK: %d, elapsed time: %g, time per point: %g. ms\n', size(kh, 1), size(kh, 2), chron, chron/length(kh(:))*1e3);
% 
% % plot wave numbers:
% ffh = fh.*ones(size(kh));
% figure(3), scatter(real(kh(:))/h/1e3, ffh(:)/h/1e6, 8, abs(imag(kh(:))), 'filled'), 
% caxis([0, 0.12]), xlim([0, 12]), ylim([0, fh(end)/h/1e6])
% xlabel('k in rad/mm'), ylabel('f in MHz')
% 
% % plot phase vel:
% ffh = fh.*ones(size(kh));
% figure(4), scatter(ffh(:)/h/1e6, 2*pi*ffh(:)./real(kh(:))/1e3, 8, abs(imag(kh(:))), 'filled'), 
% caxis([0, 0.12]), xlim([0, fh(end)/h/1e6]), ylim([0, 15]),
% ylabel('cp in mm/us'), xlabel('f in MHz')

