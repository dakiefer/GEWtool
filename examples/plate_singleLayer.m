%% compute guided ultrasonic waves in a plate
% 
% 2022 Daniel Kiefer
% Institut Langevin, Paris, France
% 

clear all
%% specify parameters:
rho = 7900; lbd = 1.1538e11; mu = 7.6923e10; % steel material
h = 1e-3; % thickness 
N = 15; % number of discretization points
k = linspace(1e-2, 12, 200)/h; % wavenumbers to solve for

%% material and waveguide description:
II = eye(3).*shiftdim(eye(3), -2); % 4th order "unit tensor"
c = lbd*II + mu*(permute(II, [1 3 4 2]) + permute(II, [1 3 2 4])); % stiffness tensor
steel.c = c; steel.rho = rho; % material 
plate = Plate(steel, h, N); % create waveguide description 

%% compute and plot:
guw = plate.Lamb; % assembles matrices for the specified waves
ff = computeW(guw, k)/2/pi; kk = k.*ones(size(ff));
figure, plot(kk(:), ff(:), '.'); ylim([0, 6e3]/h);
xlabel('wavenumber k in rad/m'), ylabel('frequency f in Hz')
