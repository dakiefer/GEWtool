% Run using: runtests()
% Compare waves in a piezoelectric plate with solutions by Kuznetsova [1]. 
%
% [1] I. E. Kuznetsova, B. D. Zaitsev, I. A. Borodina, A. A. Teplyh, V. V.
% Shurygin, and S. G. Joshi, “Investigation of acoustic waves of higher order
% propagating in plates of lithium niobate,” Ultrasonics, vol. 42, no. 1, pp.
% 179–182, Apr. 2004, doi: 10.1016/j.ultras.2004.01.006.
%
% see also: 
% https://fr.mathworks.com/help/matlab/matlab_prog/write-script-based-unit-tests.html
%
% 2024 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France

h = 1e-3; 
k = linspace(0.15, 15, 60).'/h; % wavenumber list where to solve 
N = 15;
mat0 = MaterialPiezoelectric('lithium_niobate');

%% Y-cut, propagation in X-direction  
matYX = mat0.rotateEuler(-pi/2,'x'); % passive intrinsic rotation
plate = Plate(matYX, h, N);
gew = plate.fullyCoupled; 
dat = computeW(gew, k, 10);
pic=imread('data/Kuznetsova_yCutXdir.png');   % load reference
fig = figure; clf; hold on
image([0,9],[10,0],pic); axis xy
xlabel('$\omega h / 2\pi$ in mm/us'); ylabel('$c_p$ in mm/us');
title('reference: Kuznetsova et al. 2004','Interpreter','none')
hh = plot(dat.w*h/2/pi/1e3, dat.w./dat.k/1e3, 'r.', 'MarkerSize',10,'DisplayName', 'GEWtool');
xlim([0 9]); ylim([0, 10]); 
legend(legendUnq, 'Location','southeast');

% request user to evaluate test:
assert( userTestConfirmation(fig) )

%% X-cut, propagation in Y+30°-direction 
matXY30 = mat0.rotateEuler(pi/2, 'y', (90+30)/180*pi, 'z'); % passive, intrinsic per default
plate = Plate(matXY30, h, N);
gew = plate.fullyCoupled; 
dat = computeW(gew, k, 10);
pic=imread('data/Kuznetsova_xCutY30dir.png');   % load reference
fig = figure; clf; hold on
image([1,9],[16,4],pic); axis xy
xlabel('$\omega h / 2\pi$ in mm/us'); ylabel('$c_p$ in mm/us');
title('reference: Kuznetsova et al. 2004','Interpreter','none')
hh = plot(dat.w*h/2/pi/1e3, dat.w./dat.k/1e3, 'r.', 'MarkerSize',10,'DisplayName', 'GEWtool');
xlim([1 9]); ylim([4, 16]); 
legend(legendUnq, 'Location','northeast');

% request user to evaluate test:
assert( userTestConfirmation(fig) )
