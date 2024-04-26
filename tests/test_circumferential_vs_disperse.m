% Compare circumferentially propagating waves on a hollow tube to the solutions
% computed using disperse.
% 
% Run using: runtests(). see also: 
% https://fr.mathworks.com/help/matlab/matlab_prog/write-script-based-unit-tests.html
%
% 2024 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France

rho = 2680; 
[lbd, mu] = MaterialIsotropic.wavespeed2lame(6243, 3145, rho);
mat = MaterialIsotropic('aluminum', lbd, mu, rho);
h = 1e-3; ri = 5e-3; ro = ri+h; rm = (ri+ro)/2; % inner and outer radii and more
N = 15; % discretization
cyl = CylinderCircumferential(mat, [ri,ro], N); % circumferential geometry
gew = cyl.fullyCoupled; % polarized in x-r-phi
load('data/disperse_circumferential.mat', 'disperse'); % load reference
km = linspace(1e-4, 25e3, 70); % (in rad/m) prescribe wavenumbers at middle radius rm
k = km*rm/ro; % wavenumbers at outer radius ro
dat = computeW(gew, k, 18);  % compute

%% compare against disperse
fig = figure(1); clf; hold on
hdisp = plot(2*pi*disperse.k, disperse.f, '-', 'Color', 0.7*[1,1,1], 'DisplayName','disperse');
ylim([0, 10]); xlim([0, 25]); 

% plot GEWtool result
ph = plot(dat.k*ro/rm/1e3, dat.w/2/pi/1e6, 'k.', 'MarkerSize', 8);
xlabel('wavenumber k in rad/mm'), ylabel('frequency w/2pi in MHz')
legend([ph(1), hdisp(1)], {'GEWtool', 'Disperse'}, 'Location','southeast')
title('circumferential waves: k at r = (ri + ro)/2')

% request user to evaluate test:
assert( userTestConfirmation(fig) )
