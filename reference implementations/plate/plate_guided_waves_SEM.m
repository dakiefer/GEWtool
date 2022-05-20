%% Dispersion calculation in a generally anisotropic elastic plate
% 
% 2022 - Daniel A. Kiefer
% Institut Langevin, Paris, France
% 

h = 1e-3;   % thickness in m
N = 10;     % discretization: polynomial order of interpolants
rho = 7900; lbd = 1.1538e11; mu = 7.6923e10; % steel material
mat = MaterialIsotropic('aluminum');

tic
guw = Lamb_matrices_SEM(mat, h, N);
L2 = guw.op.L2; L1 = guw.op.L1; L0 = guw.op.L0; M = guw.op.M;
fh0 = guw.np.fh0;
chron = toc; fprintf('assembling time: %g\n', chron);

%% solve for frequency and plot:
kh = linspace(1e-2, 15, 300); % wavenumber*thickness 
whn = nan(size(M, 2), length(kh)); tic 
for ii = 1:length(kh)
    kh0 = kh(ii);
    [wh2] = eig((1i*kh0)^2*L2 + (1i*kh0)*L1 + L0, -M, 'chol'); 
    whn(:,ii) = sqrt(wh2);
end
fh = whn/2/pi*fh0; fh(fh == 0) = nan;
chron = toc; fprintf('nF: %d, nK: %d, elapsed time: %g, time per point: %g. ms\n', size(fh, 1), size(fh, 2), chron, chron/length(fh(:))*1e3);

% plot wavenumbers:
kkh = kh.*ones(size(fh));
figure(1), hold on, plot(kkh(:)/h/1e3, real(fh(:))/h/1e6, '.');
xlim([0, 12]), ylim([0, 6]),
xlabel('k in rad/mm'), ylabel('f in MHz'),

ind = find(imag(fh)~=0);
figure(1), hold on, plot(kkh(ind)/h/1e3, real(fh(ind))/h/1e6, 'x');
drawnow

%% solve for wavenumbers:
fh = linspace(1e-2, 6, 300).'*1e6*h; % frequency*thickness
kh = nan(length(fh), size(M, 2)*2); tic
for ii = 1:length(fh)
    whn = 2*pi*fh(ii)/fh0; % current frequency-thickness (normalized)
    [un, khi] = polyeig(L0 + whn^2*M, L1, L2); 
    kh(ii, :) = -1i*khi;
end
chron = toc; fprintf('nF: %d, nK: %d, elapsed time: %g, time per point: %g. ms\n', size(kh, 1), size(kh, 2), chron, chron/length(kh(:))*1e3);

% % plot wave numbers:
ffh = fh.*ones(size(kh));
figure(2), hold on, scatter(real(kh(:))/h/1e3, ffh(:)/h/1e6, 8, abs(imag(kh(:))), 'filled'), 
caxis([0, 0.12]), xlim([0, 12]), ylim([0, fh(end)/h/1e6])
xlabel('k in rad/mm'), ylabel('f in MHz')

