clear all

r = [5e-3, 6e-3];
N = [20];
Nudof = 3;
mats = Material('zircaloy', 99.3e9, 0.37, 6560, 'Enu');
n = 1;
udofs = 1:3;

%% derived parameters and normalize parameters:
c0 = mats(1).tensor(1,2,1,2); h0 =(r(end)-r(1))/length(mats); % normalization parameters
rho0 = mats(1).rho; f0 = sqrt(c0/rho0)/h0; % normalization parameters

geom = Geometry(r, N, Nudof);
l1 = LayerCylindrical(mats(1), r(1), r(2), N(1));

[M, L0, L1, L2] = assembleLayers(geom, [l1], n, c0, h0, rho0);

% incorporate BCs:
[M, L0, L1, L2] = freeBCs(M, L0, L1, L2, geom, [l1], n, c0, h0);

% hl = l1.h; c0l = l1.mat.tensor(1,2,1,2); rho0l = l1.mat.rho;
% [B0, B1] = l1.tractionOp(udofs, n); B0 = B0*c0l/c0*(h0/hl)^2; B1 = B1*c0l/c0*h0/hl;
% dofBC = [1, N, N+1, 2*N, 2*N+1, 3*N]; 
% L0(dofBC, :) = B0; L1(dofBC, :) = B1; L2(dofBC, :) = 0; M(dofBC, :) = 0;

%% solve for frequency:
kh = linspace(1e-2, 10, 500); % wavenumber-thickness (solve for frequency)
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
% kkh = kh.*ones(size(fh));
% figure, plot(kkh(:), fh(:), '.');
% ylim([0, 6e3]),
% ylabel('fh in m/s'), xlabel('kh in 1')

% phase vel:
% figure, hold on
kkh = kh.*ones(size(fh));
plot(fh(:), 2*pi*fh(:)./kkh(:), '.');
xlim([0, 4e3]), ylim([0, 10e3])
xlabel('fh in m/s'), ylabel('cp in m/s')
