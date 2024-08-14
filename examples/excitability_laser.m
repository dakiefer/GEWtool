%% Laser-ultrasonic excitability of guided waves
% Plot the laser-ultrasonic excitability/detectability of guided waves in a
% bilayered plate.
% 
% 2022 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France

% % specify parameters:
mat1 = Material('steel'); mat2 = Material('chromium'); % load material data
h1 = 1e-3; h2=0.2e-3; % thickneses
N1 = 12; N2 = 8; % number of discretization points
k = linspace(1, 12e3, 200); % wavenumbers to solve for
nModes = 8;
plate = Plate({mat1, mat2}, [h1, h2], [N1, N2]); % create waveguide description 
% % compute frequencies and plot:
gew = plate.Lamb; % assembles matrices for the specified waves
dat = computeW(gew, k, nModes); 

%% plot laser-ultrasonic excitability on the "top" side, i.e., the last layer:
% Note: The excitability is normalized such that 1 corresponds to 100x the median 
% excitability. It can, hence, vary somewhat from problem to problem. To mantain
% the same normalization with different discretiations N1, N2, it is important
% to either (1) specify the number of modes to take into account as is done
% above, or (2) restrict to the solutions in a given frequency range. Otherwise,
% more and more modes are considered with increasing discretization order. 
exc = excitabilityLUS(gew, dat, 'top');
exc = 20*log10(exc); % in decibel
[exc, ind] = sort(exc(:)); % plot high excitability last (on top)
ks = dat.k(ind); ws = dat.w(ind);
figure, hold on, ylim([0, 6e6]); 
scatter(ks, ws/2/pi, 8, exc, 'filled')
colormap(flipud(colormap));
cb = colorbar; caxis([-60, 0]);
cb.Label.Interpreter = 'latex';  cb.Label.String = '$u_x u_r$ at outer surface in dB';
xlabel('wavenumber k in rad/m'), ylabel('frequency f in Hz')
title(sprintf('zircaloy plate'))
