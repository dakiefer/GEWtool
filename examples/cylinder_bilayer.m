%% compute axially guided ultrasonic waves in a cylinder
% waves in a cylinder with two layeres of different material
% 
% 2022 Daniel Kiefer
% Institut Langevin, Paris, France
% 

% specify parameters:
steel = Material('steel');
zirc = Material('zircaloy', 99.3e9, 0.37, 6560, 'Enu'); 
r = [5e-3, 5.37e-3 6e-3]; % radial coordinates of layer interfaces
N = [15, 16]; % number of discretization points
n = 0; % circumferential order
k = linspace(1e-2, 15, 200)/(r(end)-r(1)); % wavenumber-thickness (solve for frequency)

% compute
cyl = Cylinder([steel, zirc], r, N); % create waveguide description 
guw = cyl.fullyCoupled(n); % waves of circumferential order n
ff = computeW(guw, k)/2/pi; kk = k.*ones(size(ff));

% plot
figure, plot(kk(:), ff(:), '.'); ylim([0, 6e3]/(r(end)-r(1)));
xlabel('wavenumber k in rad/m'), ylabel('frequency f in Hz')
