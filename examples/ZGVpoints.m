%% Compute ZGV points of guided elastic waves in a plate
% Three different computational techniques are implemented in GEWtool and
% these are showcased in the following.
%
% For details refer to:
% D. A. Kiefer, B. Plestenjak, H. Gravenkamp, and C. Prada, "Computing
% zero-group-velocity points in anisotropic elastic waveguides: Globally and
% locally convergent methods," The Journal of the Acoustical Society of America,
% vol. 153, no. 2, pp. 1386–1398, Feb. 2023, doi: 10.1121/10.0017252.
%
% 2022 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France

% % specify parameters:
mat = Material('steel_austenitic'); % orthotropic 
mat = mat.rotateEuler(180/180*pi, 90/180*pi, 0); % this direction is fun
h = 1e-3; % thickness 
N = 25; % number of discretization points
wmax = 2*pi*15e6; % maximum frequency for plotting and ZGV-search
k = linspace(1e-2, 40, 150)/h; % wavenumbers for plotting dispersion curves
plate = Plate(mat, h, N); % create waveguide description 
gew = plate.LambA; % assembles matrices for the specified waves
dat = computeW(gew, k); % dispersion curves will be plotted for reference 

% % here are some ZGV points with 2 digits accuracy (serve as initial guess): 
w0 = [0.29    0.52    0.69    0.72    0.75]*1e8;
k0 = [0.34    0.35    0.71    0.25    0.07]*1e4;

% % plot dispersion curves for reference: 
fig = figure; clf, hold on, grid off;
plot(dat.k*h, dat.w*h/2/pi/1e3, 'Color', [.5, .5, .5], 'HandleVisibility', 'off');
ylim([0, wmax/2/pi/1e6]), xlim([0, wmax/2/pi/1e6]), legend('Location','southeast')
set(fig,'defaulttextinterpreter','latex'), set(fig, 'DefaultLegendInterpreter', 'latex')
xlabel('wavenumber-thickness $kh$ in rad'), ylabel('frequency-thickness $fh$ in MHz$\cdot$mm')
drawnow;

%% Newton-kind iteration: 
% This method is implemented in computeZGV(). It is super fast but needs initial
% guesses:
fprintf('\n\n++ Newton-type iteration: ++\n')
fprintf('Search close to provided initial guesses:\n')
tic, zgv5 = computeZGV(gew, w0, k0); toc
% print initial guesses and computed values 
disp('initial frequencies:'), disp(w0/2/pi)
disp('converged frequencies:'), disp(zgv5.w.'/2/pi)

% Instead of initial guesses w0, k0, you can also provide the dispersion data
% "dat", computeZGV() will then search where the group velocity changes sign:
fprintf('\n\n++ Newton-type iteration: ++\n')
fprintf('Search where group velocity changes sign:\n')
tic, zgvNewton = computeZGV(gew, dat); timing = toc;
nZGV = length(zgvNewton.w(zgvNewton.w < wmax)); 
fprintf('Computed %d ZGV points in %g s.\n', nZGV, timing); 

figure(fig);
plot(zgvNewton.k(:)*h, zgvNewton.w(:)*h/2/pi/1e3, 'rx', 'MarkerSize', 10, 'DisplayName', 'Newton method');
drawnow;

%% Scanning the ZGV points
% This methods is implemented in computeZGVScan(). It does not need initial
% guesses and is likely to locate all ZGV points but it is substantially slower
% than the Newton-type iteration. To speed up the computation, provide the
% maximum frequency to search at. Even more effective is to provide the paramter
% kMax (and optionally kStart), these define the wavenumber search interval
% [kStart, kMax].
fprintf('\n\n++ Scanning method: ++\n')
clear opts, opts.kEnd = 12/h; opts.show = false; % unit distance is np.h0
tic, zgvScan = computeZGVScan(gew, wmax, opts); timing = toc;
nZGV = length(zgvScan.w(zgvScan.w < wmax));
fprintf('Computed %d ZGV points in %g s.\n', nZGV, timing); 

figure(fig);
plot(zgvScan.k(:)*h, zgvScan.w(:)*h/2/pi/1e3, 'ksquare', 'DisplayName', 'Scanning method', 'MarkerSize', 8);
if isfield(zgvScan, 'k0s'), xline(zgvScan.k0s*h, 'HandleVisibility', 'off'); end
drawnow;

%% Direct method
% This method is implemented in computeZGVDirect(). It does not need initial
% guesses and guarantees to find all ZGV points as long as the matrices defining
% the problem are not too big. This is rather slow and should not be used for
% matrices bigger than about 40x40. 
fprintf('\n\n++ Direct method: ++\n')
tic, zgvDirect = computeZGVDirect(gew); timing = toc;
nZGV = length(zgvDirect.w(zgvDirect.w < wmax));
fprintf('Computed %d ZGV points in %g s.\n', nZGV, timing); 

figure(fig);
plot(zgvDirect.k(:)*h, zgvDirect.w(:)*h/2/pi/1e3, 'k.', 'MarkerSize', 8, 'DisplayName', 'Direct method');
drawnow;
