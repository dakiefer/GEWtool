% specify parameters:
clear all
r = [20e-3, 20.5e-3, 21e-3];
N = [30, 30];
Nudof = [3, 3];
mat = Material('steel');
% mat = Material('zircaloy', 99.3e9, 0.37, 6560, 'Enu'); 
mats = [mat, mat];
n = 0;
kh = linspace(1e-3, 5, 200); % wavenumber-thickness (solve for frequency)

%% one layer:
c0 = mats(1).tensor(1,2,1,2); h0 = r(3) - r(1); % normalization parameters
rho0 = mats(1).rho; fh0 = sqrt(c0/rho0); % normalization parameters
geom = Geometry([r(1), r(3)], 40, Nudof);
l1 = LayerCylindrical(mats(1), r(1), r(3), 40);
[M, L0, L1, L2] = assembleLayers(geom, l1, n, c0, h0, rho0);
[M, L0, L1, L2] = freeBCs(M, L0, L1, L2, geom, l1, n, c0, h0);

% solve and plot:
whn = solveWithKh(kh, M, L0, L1, L2);
fh = whn/2/pi*fh0; 
kkh = kh.*ones(size(fh));
figure, hold on, plot(kkh(:)/h0, fh(:)/h0, 'gx');

%% two layers of same material and same total thickness as one layer:
geom = Geometry(r, N, Nudof);
l1 = LayerCylindrical(mats(1), r(1), r(2), N(1));
l2 = LayerCylindrical(mats(2), r(2), r(3), N(2));
[M, L0, L1, L2] = assembleLayers(geom, [l1, l2], n, c0, h0, rho0);
[M, L0, L1, L2] = freeBCs(M, L0, L1, L2, geom, [l1, l2], n, c0, h0);

% solve and plot:
whn = solveWithKh(kh, M, L0, L1, L2);
fh = whn/2/pi*fh0; 
kkh = kh.*ones(size(fh));
plot(kkh(:)/h0, fh(:)/h0, 'k.');

xlim([0, 5e3]), ylim([0, 4e6])
xlabel('k in rad/m'), ylabel('f in Hz')
legend({'one', 'two'})
