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
mat = mat.rotateEuler(90/180*pi, 90/180*pi, 0); % this direction is fun
h = 1e-3; % thickness 
N = 25; % number of discretization points
fmax = 15e6; % maximum frequency for plotting and ZGV-search
k = linspace(1e-2, 40, 150)/h; % wavenumbers for plotting dispersion curves
plate = Plate(mat, h, N); % create waveguide description 
gew = plate.LambA; % assembles matrices for the specified waves
dat = computeW(gew, k); % dispersion curves will be plotted for reference 
dat.cg = groupVel(gew, dat); % compute group velocity for plotting

% % here are some ZGV points with 2 digits accuracy (serve as initial guess): 
w0 = [0.29    0.52    0.69    0.72    0.75]*1e8;
k0 = [0.34    0.35    0.71    0.25    0.07]*1e4;

%% Newton-kind iteration: 
% This method is implemented in computeZGV(). It is super fast but needs initial
% guesses:
fprintf('\n\n++ Newton-type iteration: ++\n')
fprintf('Search close to provided initial guesses:\n')
tic, zgv5 = computeZGV(gew, w0, k0); toc
w0initial = w0, wConverged = zgv5.w.' % print initial guesses and computed values

% Instead of initial guesses w0, k0, you can also provide the dispersion data
% "dat", computeZGV() will then search where the group velocity changes sign:
fprintf('\n\n++ Newton-type iteration: ++\n')
fprintf('Search where group velocity changes sign:\n')
tic, zgv = computeZGV(gew, dat); toc
nZGV = length(zgv.w(zgv.w < 2*pi*fmax)) % print number of ZGV points 
plotZGVs(dat, zgv, fmax), title('Newton-type method'), drawnow


%% Scanning the ZGV points
% This methods is implemented in computeZGVScan(). It does not need initial
% guesses and is likely to locate all ZGV points but it is substantially slower
% than the Newton-type iteration. To speed up the computation, provide the
% maximum frequency to search at. Even more effective is to provide the paramter
% kMax (and optionally kStart), these define the wavenumber search interval
% [kStart, kMax].
fprintf('\n\n++ Scanning method: ++\n')
clear opts, opts.kMax = 10/gew.np.h0; % unit distance is np.h0
tic, zgv = computeZGVScan(gew, opts); toc
nZGV = length(zgv.w(zgv.w < 2*pi*fmax)) % print number of ZGV points 
plotZGVs(dat, zgv, fmax), title('Scaning method'), drawnow

%% Direct method
% This method is implemented in computeZGVDirect(). It does not need initial
% guesses and guarantees to find all ZGV points as long as the matrices defining
% the problem are not too big. This is rather slow and should not be used for
% matrices bigger than about 40x40. 
fprintf('\n\n++ Direct method: ++\n')
tic, zgv = computeZGVDirect(gew); toc
nZGV = length(zgv.w(zgv.w < 2*pi*fmax)) % print number of ZGV points 
plotZGVs(dat, zgv, fmax), title('Direct method'), drawnow


%% function to plot dispersion curves with ZGV points
function plotZGVs(dat, zgv, fmax)
    fh = figure; pos = fh.Position; fh.Position = [pos(1), pos(2)-pos(4)*0.8, pos(3), pos(4)*1.8];

    subplot(2,1,1), hold on, ylim([0, fmax]); xlim([0, 25e3])
    plot(dat.k, dat.w/2/pi, '-'); ax = gca; ax.ColorOrderIndex = 1; 
    xlabel('wavenumber k in rad/m'), ylabel('frequency f in Hz')
    plot(zgv.k(:), zgv.w(:)/2/pi, 'r*');
    
    subplot(2,1,2), hold on, xlim([0, fmax]); 
    plot(real(dat.w)/2/pi, real(dat.cg), '-'); ax = gca; ax.ColorOrderIndex = 1; 
    xlabel('frequency f in Hz'), ylabel('group vel in m/s')
    plot(zgv.w(:)/2/pi, zeros(size(zgv.w(:))), 'r*')
    subplot(2,1,1);
end
