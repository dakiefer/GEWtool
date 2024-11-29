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
gew = plate.LambSA; tic;        % choose S+A Lamb waves (assembles matrices)
dat = computeW(gew, k, 4); toc; % solve and save 4 modes (argument optional)

figure(1); clf; hold on
plot(dat(1).k/1e3, dat(1).w/2/pi/1e6,'SeriesIndex',1,'DisplayName','symmetric');
plot(dat(2).k/1e3, dat(2).w/2/pi/1e6,'SeriesIndex',2,'DisplayName','anti-symm.');
ylim([0, dat(1).w(end,1)/2/pi/1e6*1.1]); 
xlabel('wavenumber k in rad/mm'), ylabel('frequency f in MHz')
legend(legendUnq, 'Location', 'southeast')
title(sprintf('Lamb waves in %gmm %s', h/1e-3, mat.name))

% % Energy velocity (identical to the group velocity)
% We plot the "axial" component cex of the energy velocity vector ce = [cex, cey], 
% i.e., the one that is along the wave vector k. For an isotropic plate this is
% the only nonzero component.
figure(2); clf; hold on; 
plot(dat(1).w/2/pi/1e6, dat(1).cex/1e3, 'SeriesIndex',1,'DisplayName','symmetric');
plot(dat(2).w/2/pi/1e6, dat(2).cex/1e3, 'SeriesIndex',2,'DisplayName','anti-symm.');
xlim([0, 6]); 
xlabel('frequency f in MHz'), ylabel('energy velocity ce in mm/us')
legend(legendUnq, 'Location', 'southeast')
title(sprintf('Lamb waves in %gmm %s', h/1e-3, mat.name))
