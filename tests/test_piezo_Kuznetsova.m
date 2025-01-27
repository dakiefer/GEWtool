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
matYX = mat0.rotateEuler(pi/2,'x');
% matXY30 = mat0.rotateEuler(-pi/2,'y',pi+30/180*pi,'x');
matXY30 = mat0.rotateEuler(pi/2,'y',(30+90)/180*pi,'z');

%% Y-cut, propagation in X-direction  
plate = Plate(matYX, h, N);
gew = plate.polarization(1:3, 0); 
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
plate = Plate(matXY30, h, N);
gew = plate.polarization(1:3, 0); 
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
