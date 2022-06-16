%% compute ZGV points of guided elastic waves in a plate
% 
% 2022 Daniel Kiefer
% Institut Langevin, Paris, France
% 

%% specify parameters:
mat = Material('triclinic');
h = 1e-3; % thickness 
N = 15; % number of discretization points
k = linspace(1e-2, 40, 600)/h; % wavenumbers to plot dispersion curves

%% material and waveguide description:
plate = Plate(mat, h, N); % create waveguide description 
gew = plate.fullyCoupled; % assembles matrices for the specified waves
dat = computeW(gew, k); % dispersion curves will be plotted as reference 

%% compute ZGV points
cg = groupVel(gew, dat);
sigChange = diff(sign(real(cg)),1,2); % detect where the sign changes
w0 = dat.w(find(sigChange));
k0 = dat.k(find(sigChange));
addpath('~/Projekte/ZGVProjekt/ZGVcomputation/code/')
tic
zgv = computeZGVCloseTo(gew, w0, k0);
% zgv = computeZGV(gew);
toc

% plot
figure(1), clf, hold on, ylim([0, 10e3]/h);
plot(dat.k.', dat.w.'/2/pi, 'k'); ax = gca; ax.ColorOrderIndex = 1; 
xlabel('wavenumber k in rad/m'), ylabel('frequency f in Hz')
plot(zgv.k(:), zgv.w(:)/2/pi, 'r*');

figure(2), clf, hold on, xlim([0, 10e3]/h); 
plot(real(dat.w).'/2/pi, real(cg).', 'k'); ax = gca; ax.ColorOrderIndex = 1; 
xlabel('frequency f in Hz'), ylabel('group vel in m/s')
plot(zgv.w(:)/2/pi, zeros(size(zgv.w(:))), 'r*')


% zgv1 = computeZGVNewton(gew, 2.7e6*2*pi, 2000);
% plot(zgv1.k(:), zgv1.w(:)/2/pi, 'r*')
% zgv1 = computeZGVNewton(gew, 4.6e6*2*pi, 1300);
% plot(zgv1.k(:), zgv1.w(:)/2/pi, 'r*')

