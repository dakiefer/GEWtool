%% Compute Lamb waves in a plate
% This example shows how to compute Lamb waves. The symmetric and anti-symmetric
% waves are computed separately. This is faster and the modes are correctly
% orderd.
%
%  _↑z____________________________________ top surface
%   →x     --> k  guided wave      c,rho   plate (stiffness c, mass density rho)
%  _______________________________________ bottom surface
%
% 2024-2025 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France

mat = MaterialIsotropic('steel'); % load from database (or create your own)
h = 1e-3;                        % thickness in m
N = 12;                          % number of nodes (dictates accuracy)
k = linspace(1e-2, 15, 100)/h;   % wavenumbers to solve for
plate = Plate(mat, h, N);        % create waveguide description 
gew = plate.LambSA; tic;         % choose S+A Lamb waves (assembles matrices)
dat = computeW(gew, k, 6); toc;  % solve and save 6 modes (argument optional)

% % plot frequency-wavenumber dispersion curves:
figure(1); clf; 
plot(dat);
title(sprintf('%s, %g mm thick', mat.name, h/1e-3))

% % Energy velocity (identical to the group velocity)
% We plot the axial component cex = ce(:,:,1) of the energy velocity vectors ce, 
% i.e., the one that is along the wave vector k. For an isotropic plate this is
% the only non-zero component.
% The property dat.ce is an alias to either (i) energyVelVec(dat) if the plate
% is anisotropic or (ii) energyVelAxial(dat) if the plate is isotropic. These
% functions can be found in the GEWdat directory. You can explore this directory
% to see what postprocessing functions are available.
figure(2); clf; hold on; 
plot(dat(1).w/2/pi/1e6, dat(1).ce/1e3, 'SeriesIndex',1,'DisplayName','symmetric');
plot(dat(2).w/2/pi/1e6, dat(2).ce/1e3, 'SeriesIndex',2,'DisplayName','anti-symm.');
xlim([0, 1]*dat(1).w(end,1)/2/pi/1e6); 
xlabel('frequency f in MHz'), ylabel('energy velocity ce in mm/us')
legend(legendUnq, 'Location', 'southeast')
title(sprintf('%s, %g mm thick', mat.name, h/1e-3))
