%% Compute circumferentially guided waves in a tube/rod
% Solve for waves travelling circumferentially around a circular tube or rod.
% The waves might segregate into two families: ur-uphi-polarized (Lamb-like) and 
% ux-polarized (SH-like).
% 
% Literature: 
% J. Qu, Y. Berthelot, and Z. Li, “Dispersion of Guided Circumferential Waves in
% a Circular Annulus,” in Review of Progress in Quantitative Nondestructive
% Evaluation: Volume 15A, D. O. Thompson and D. E. Chimenti, Eds., Boston, MA:
% Springer US, 1996, pp. 169–176. doi: 10.1007/978-1-4613-0383-1_21.
% 
% 2024 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France

mat = MaterialIsotropic('aluminum');
h = 1e-3; ri = 1e-3; ro = ri+h; % thickness, inner and outer radii
N = 10; % number of nodes on radial coordinate
cyl = CylinderCircumferential(mat, [ri,ro], N); % cylindrical geometry + circumf. prop.
k = linspace(1e-4, 20e3, 150); % (in rad/m) prescribe wavenumbers at outer radius ro

% compute Lamb-like
gew = cyl.Lamb;                 % polarized in r-phi (radial-angular directions)
datL = computeW(gew, k, 12);    % compute
% compute SH-like
gew = cyl.sh;                   % polarized in x (axial direction)
datSH = computeW(gew, k, 8);    % compute

% plot frequency-wavenumber dispersion 
fig = figure(1); clf; hold on;
phSH = plot(datSH.k/1e3, datSH.w/2/pi/1e6, 'Color', "#5EC962");
phL = plot(datL.k/1e3, datL.w/2/pi/1e6, 'Color', "#3B518B");
ylim([0, 6]); xlim([0, 10]);
xlabel('wavenumber ko in rad/mm'), ylabel('frequency w/2pi in MHz')
legend([phL(1), phSH(1)], {'r-phi-polarized', 'x-polarized'}, 'Location','southeast')
title(sprintf('circumferential waves (ri: %g mm, ro: %g mm)', ri/1e-3, ro/1e-3));
