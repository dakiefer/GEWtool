% Compare axial waves in a solid rod to the solutions computed using Disperse.
%
% Run using: runtests(). see also: 
% https://fr.mathworks.com/help/matlab/matlab_prog/write-script-based-unit-tests.html
%
% 2024 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France

% rod parameters and other definitions
rho = 2680; ct = 3145; cl = 6243; % aluminum material parameters 
[lbd, mu] = MaterialIsotropic.wavespeed2lame(cl, ct, rho);
mat = MaterialIsotropic('aluminum', lbd, mu, rho);
r = 1e-3;                            % radius r
N = 15;                              % number of nodes (dictates accuracy)
k = linspace(1e-4, 15, 100)/r;       % wavenumbers to solve for
cyl = Cylinder(mat, [0, r], N);      % create waveguide description 
load("data/disperse_rod_axial.mat","disperse"); % load reference

%% compare against disperse
% compute longitudinal and flexural modes:
long = cyl.longitudinal;  % assembles matrices 
datL = computeW(long, k);
flex = cyl.flexural(1); % first order flexural modes
datF = computeW(flex, k);

% plot frequency-wavenumber dispersion against reference:
fig = figure(1); clf; hold on
rdisp = 1e-3;
hdisp = plot(1e3*disperse.k*2*pi*rdisp, 2*pi*1e6*disperse.f*rdisp/ct, '-', 'Color',0.7*[1,1,1]);
hL = plot(datL.k*r, datL.w*r/ct, 'k.', 'MarkerSize', 8); % plot longitudinal in blue 'Color', "#3B518B"
hF = plot(datF.k*r, datF.w*r/ct, 'k.', 'MarkerSize', 8); % plot torsional in red 'Color', "#5EC962"
ylim([0, 10]); xlim([0, 10]);
xlabel('$k r$'), ylabel('$\omega r / c_t$')
legend([hL(1), hdisp(1)], {'GEWtool','disperse'}, 'Location','south east')
title(sprintf('%s rod', mat.name))

% request user to evaluate test:
assert( userTestConfirmation(fig) )
