%% Compute Lamb waves in a plate
% This example shows how to compute Lamb waves. The symmetric and anti-symmetric
% waves are computed separately. This is faster and the modes are correctly
% orderd.
% 
% 2024 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France

mat = MaterialIsotropic('steel');   % load from database (or create your own)
h = 1e-3;                        % thickness in m
N = 12;                          % number of nodes (dictates accuracy)
k = linspace(1e-2, 12, 100)/h;   % wavenumbers to solve for
plate = Plate(mat, h, N);        % create waveguide description 
gews = plate.LambSA; tic;        % choose S+A Lamb waves (assembles matrices)
dat = computeW(gews, k, 4); toc; % solve and save 4 modes (argument optional)

figure(1); clf; hold on
hS = plot(dat(1).k/1e3, dat(1).w/2/pi/1e6, 'Color', "#3B518B"); % plot symmetric
hA = plot(dat(2).k/1e3, dat(2).w/2/pi/1e6, 'Color', "#5EC962"); % plot anti-symmetric
ylim([0, 6]); 
xlabel('wavenumber k in rad/mm'), ylabel('frequency f in MHz')
legend([hS(1), hA(1)], {'symmetric', 'anti-symmetric'}, 'Location', 'southeast')
title(sprintf('Lamb waves in %gmm %s', h/1e-3, mat.name))
