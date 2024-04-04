%% Compute waves in a circular rod
% Solve for waves travelling axially along a non-hollow circular rod. GEWtool
% automatically switches the Finite Element integration scheme to avoid the
% singularity at r=0 that appears in the 1/r-term. Only waves of 0th 
% circumferential order are computed. These segregate into two families of
% waves: torsional (uPhi-polarized) and longitudinal (ux-ur-polarized). 
% 
% 2024 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France

mat = MaterialIsotropic('steel');    % load from material/database
r = 1e-3;                            % radius r
N = 10;                              % number of nodes (dictates accuracy)
k = linspace(1e-4, 15, 200)/r;       % wavenumbers to solve for
cyl = Cylinder(mat, [0, r], N); tic; % create waveguide description 

% compute torsional waves (uPhi displacements only)
tors = cyl.torsional;     % assembles matrices
datT = computeW(tors, k);

% compute longitudinal waves (ux and ur displacements only)
long = cyl.longitudinal;  % assembles matrices 
datL = computeW(long, k); toc

% plot frequency-wavenumber dispersion 
figure(1); clf; hold on
hT = plot(datT.k/1e3, datT.w/2/pi/1e6, 'Color', "#5EC962"); % plot torsional in red
hL = plot(datL.k/1e3, datL.w/2/pi/1e6, 'Color', "#3B518B"); % plot longitudinal in blue
ylim([0, 6]); % frequency range
xlabel('wavenumber k in rad/mm'), ylabel('frequency f in MHz')
legend([hT(1), hL(2)], {'torsional', 'longitudinal'}, 'Location','south east')
title(sprintf('rod of radius %g mm', r/1e-3))

% plot phase velocity-frequency dispersion 
figure(2); clf; hold on; 
hT = plot(datT.w/2/pi/1e6, datT.w./datT.k/1e3, 'Color', "#5EC962"); % plot torsional in red
hL = plot(datL.w/2/pi/1e6, datL.w./datL.k/1e3, 'Color', "#3B518B"); % plot longitudinal in blue
xlim([0, 6]); ylim([0, mat.cl*2/1e3]) 
ylabel('phase velocity in mm/us'), xlabel('frequency f in MHz')
legend([hT(1), hL(2)], {'torsional', 'longitudinal'}, 'Location','south east')
title(sprintf('rod of radius %g mm', r/1e-3))
