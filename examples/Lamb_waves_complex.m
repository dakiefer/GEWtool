%% Compute Lamb waves in a plate in complex wavenumbers
% Compute the complex wavenumber spectrum of Lamb waves in a plate. The
% symmetric and anti-symmetric waves are computed separately. This is faster and
% the modes are correctly orderd.
% 
% 2024 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France

mat = MaterialIsotropic('steel'); % load from database (or create your own)
h = 1e-3;                         % thickness in m
N = 15;                           % number of nodes (dictates accuracy)
w = 2*pi*linspace(1e-2, 7, 1000).'*1e6; % frequencies where to compute wavenumbers k
plate = Plate(mat, h, N);         % create waveguide description 
gews = plate.LambSA; tic;         % choose S+A Lamb waves (assembles matrices)
% linearizeInK(gews);             % optional: this makes the computation faster
dat = computeK(gews, w); toc;     % solve 

figure(1); clf; hold on
hS = plot3(real(dat(1).k(:))/1e3, imag(dat(1).k(:))/1e3, dat(1).w(:)/2/pi/1e6, '.', 'Color', "#0072BD"); 
hA = plot3(real(dat(2).k(:))/1e3, imag(dat(2).k(:))/1e3, dat(2).w(:)/2/pi/1e6, '.', 'Color', "#D95319"); 
xlim([0, 12]), ylim([-10.5, 10.5]), view(22, 18)
xlabel('Re(k) in rad/mm'), ylabel('Im(k) in rad/mm'), zlabel('f in MHz')
legend([hS(1), hA(1)], {'symmetric', 'anti-symmetric'}, 'Location', 'southeast')
title(sprintf('Lamb waves in %gmm %s', h/1e-3, mat.name))
