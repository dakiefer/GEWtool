%% Compute ZGV points of guided elastic waves in a plate
% 
% 2022 Daniel A. Kiefer
% Institut Langevin, ESPCI Paris, France
% 

% % specify parameters:
mat = Material('steel_austenitic');
mat = mat.rotateEuler(90/180*pi, 90/180*pi, 0);
h = 1e-3; % thickness 
N = 40; % number of discretization points
k = linspace(1e-2, 25, 150)/h; % wavenumbers to plot dispersion curves

% % material and waveguide description:
plate = Plate(mat, h, N); % create waveguide description 
gew = plate.LambA; % assembles matrices for the specified waves
dat = computeW(gew, k); % dispersion curves will be plotted as reference 
dat.cg = groupVel(gew, dat); % compute group velocity for plotting

% w0 = [0.2883    0.5216    0.6940    0.7191    0.7512]*1e8;
% k0 = [0.3392    0.3496    0.7085    0.2491    0.0680]*1e4;
%% compute ZGV points
tic
% clear opts, opts.kMax = 20/gew.np.h0; opts.show = false;
fmax = 25e6;
zgv = computeZGV(gew, dat);
toc

nZGV = length(zgv.w(zgv.w < 2*pi*fmax))
% plot
figure(1), clf, hold on, ylim([0, fmax]); xlim([0, 25e3])
plot(dat.k.', dat.w.'/2/pi, '-k'); ax = gca; ax.ColorOrderIndex = 1; 
xlabel('wavenumber k in rad/m'), ylabel('frequency f in Hz')
plot(zgv.k(:), zgv.w(:)/2/pi, 'r*');

% plot wave speeds:
wmax = 2*pi*fmax;
waveSpeeds = vertcat(gew.lay.mat.wavespeeds); 
kMax = wmax./waveSpeeds.';
hold on
for i=1:length(kMax)
    plot([0;kMax(i)], [0;wmax/2/pi], 'Color', [.7 .7 .7], 'LineWidth',.4);
end
% 
% figure(2), clf, hold on, xlim([0, 10e3]/h); 
% plot(real(dat.w).'/2/pi, real(dat.cg).', 'k'); ax = gca; ax.ColorOrderIndex = 1; 
% xlabel('frequency f in Hz'), ylabel('group vel in m/s')
% plot(zgv.w(:)/2/pi, zeros(size(zgv.w(:))), 'r*')
