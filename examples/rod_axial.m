%% Compute waves in a circular rod
% Solve for waves travelling axially along a solid circular rod. GEWtool
% automatically switches the Finite Element integration scheme to avoid the
% singularity at r=0 that appears in the 1/r-term. The waves of 0th 
% circumferential order segregate into two families of
% waves: torsional (uPhi-polarized) and longitudinal (ux-ur-polarized). 
% 
% Literature: 
% [1] H. Gravenkamp, “Numerical methods for the simulation of ultrasonic guided
% waves,” Doctoral Dissertation, BAM Bundesanstalt für Materialforschung und
% -prüfung, Berlin, 2014.
% [2] J. Zemanek Jr., “An Experimental and Theoretical Investigation of Elastic
% Wave Propagation in a Cylinder,” The Journal of the Acoustical Society of
% America, vol. 51, no. 1B, pp. 265–283, Jan. 1972, doi: 10.1121/1.1912838.
% 
% 2024 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France

mat = MaterialIsotropic('steel');    % load from material/database
r = 1e-3;                            % radius r
N = 10;                              % number of nodes (dictates accuracy)
k = linspace(1e-4, 15, 200)/r;       % wavenumbers to solve for
cyl = Cylinder(mat, [0, r], N); tic; % create waveguide description 

% compute longitudinal waves (ux and ur displacements only)
long = cyl.longitudinal; tic;  % assembles matrices 
datL = computeW(long, k); toc,
% compute torsional waves (uPhi displacements only)
tors = cyl.torsional; tic;     % assembles matrices
datT = computeW(tors, k); toc
% compute first-order flexural waves (all displacements are coupled)
flex = cyl.flexural(1); tic;   % circumferential order n = 1
datF = computeW(flex, k); toc

% plot frequency-wavenumber dispersion 
figure(1); clf; hold on
plot(datT.k/1e3, datT.w/2/pi/1e6,'SeriesIndex',3,'DisplayName','torsional');
plot(datF.k/1e3, datF.w/2/pi/1e6,'SeriesIndex',2,'DisplayName','1st flexural');
plot(datL.k/1e3, datL.w/2/pi/1e6,'SeriesIndex',1,'DisplayName','longitudinal');
ylim([0, 5]); % frequency range
xlabel('wavenumber k in rad/mm'), ylabel('frequency f in MHz')
legend(legendUnq, 'Location','south east')
title(sprintf('rod of radius %g mm', r/1e-3))

% % plot phase velocity-frequency dispersion 
figure(2); clf; hold on; 
plot(datT.w/2/pi/1e6, datT.cp/1e3,'SeriesIndex',3,'DisplayName','torsional');
plot(datF.w/2/pi/1e6, datF.cp/1e3,'SeriesIndex',2,'DisplayName','1st flexural');
plot(datL.w/2/pi/1e6, datL.cp/1e3,'SeriesIndex',1,'DisplayName','longitudinal');
xlim([0, 5]); ylim([0, mat.cl*2/1e3]) % axes limits
ylabel('phase velocity in mm/us'), xlabel('frequency f in MHz')
legend(legendUnq, 'Location','south east')
title(sprintf('rod of radius %g mm', r/1e-3))
