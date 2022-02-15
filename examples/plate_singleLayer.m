%% compute guided ultrasonic waves in a plate
% 
% 2022 Daniel Kiefer
% Institut Langevin, Paris, France
% 

% specify parameters:
steel = Material('steel'); % material
h = 1e-3; % thickness 
N = 15; % number of discretization points
k = linspace(1e-2, 12, 200)/h; % wavenumbers to solve for
plate = Plate(steel, h, N); % create waveguide description 

% compute and plot
figure, hold on; ax = gca;
xlabel('wavenumber k in rad/m'), ylabel('frequency f in Hz')
guw = plate.Lamb(0); % coupled Lamb and SH-polarized waves
ff = computeW(guw, k)/2/pi; kk = k.*ones(size(ff));
plot(kk(:), ff(:), '.'); ylim([0, 6e3]/h);
