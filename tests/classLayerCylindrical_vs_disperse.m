clear all
load('../reference implementations/cylindrical/comparison_disperse/disperse.mat')
figure, hold on
plot(disperse.k/1e-3, disperse.f*1e6, 'g.')
xlim([0, 0.5]), ylim([0, 200])

a = 20; b = 21;
rho=7932; mu=rho*3260^2; lbd=rho*5960^2-2*mu;	
mat = Material('noname', lbd, mu, rho, 'lame');
N = 25;
Nudof = 3;
udofs = 1:3;

%% problem definition
c0 = mat.c(1,2,1,2); h0 = b - a; % normalization parameters
rho0 = mat.rho; f0 = sqrt(c0/rho0)/h0; % normalization parameters
geom = Geometry([a, b], N, Nudof);
l1 = LayerCylindrical(mat, a, b, N);

%% assemble problem
n = 0;
[M, L0, L1, L2] = assembleLayers(geom, l1, n, c0, h0, rho0);
[M, L0, L1, L2] = freeBCs(M, L0, L1, L2, geom, l1, n, c0, h0);

%% solve for frequency:
kh = linspace(1e-3, 0.5, 100); % wavenumber-thickness (solve for frequency)
tic 
whn = nan(size(M, 2), length(kh));
for ii = 1:length(kh)
    [wh2] = polyeig((1i*kh(ii))^2*L2 + (1i*kh(ii))*L1 + L0, M); % does not work properly with eig()
    whn(:,ii) = sqrt(wh2);
end
fh = real(whn/2/pi*f0*h0); % fh(fh == 0) = nan;
chron = toc;
fprintf('nF: %d, nK: %d, elapsed time: %g, time per point: %g. ms\n', size(fh, 2), size(fh, 1), chron, chron/length(fh(:))*1e3);

% wavenumbers:
kkh = kh.*ones(size(fh));
hold on, plot(kkh(:)/h0, fh(:)/h0, '.');
xlim([0, 0.5]), ylim([0, 200])
ylabel('f in Hz'), xlabel('k in rad/m')
legend({'disperse', 'n = 0'})

%% assemble problem for n = 1
n = 1;
[M, L0, L1, L2] = assembleLayers(geom, l1, n, c0, h0, rho0);
[M, L0, L1, L2] = freeBCs(M, L0, L1, L2, geom, l1, n, c0, h0);

%% solve for frequency:
kh = linspace(1e-3, 0.5, 100); % wavenumber-thickness (solve for frequency)
tic 
whn = nan(size(M, 2), length(kh));
for ii = 1:length(kh)
    [wh2] = polyeig((1i*kh(ii))^2*L2 + (1i*kh(ii))*L1 + L0, M); % does not work properly with eig()
    whn(:,ii) = sqrt(wh2);
end
fh = real(whn/2/pi*f0*h0); % fh(fh == 0) = nan;
chron = toc;
fprintf('nF: %d, nK: %d, elapsed time: %g, time per point: %g. ms\n', size(fh, 2), size(fh, 1), chron, chron/length(fh(:))*1e3);

% wavenumbers:
kkh = kh.*ones(size(fh));
hold on, plot(kkh(:)/h0, fh(:)/h0, 'k.');
xlim([0, 0.5]), ylim([0, 200])
ylabel('f in Hz'), xlabel('k in rad/m')
legend({'disperse', 'n = 0', 'n = 1'})
