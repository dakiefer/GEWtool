%% compute axially guided ultrasonic waves in bilayered cylinder
% Waves in a cylinder with two layeres of different material.
% Two different deposited layer thicknesses are compared.
% 
% 2022 Daniel Kiefer
% Institut Langevin, Paris, France
% 

% specify parameters:
ri = 5e-3; ro = 5.57e-3; % inner and outer radii
h1 = 10e-6; h2 = 30e-6; % two different deposited layer thicknesses 
N = [20, 10]; % number of discretization points
n = 0; % circumferential order
k = linspace(1e-3, 30, 400)/(r(end)-r(1)); % wavenumber-thickness (solve for frequency)
zirc = Material('zircaloy_aniso');
chrom = Material('chromium');

% first deposited layer thickness
cyl = Cylinder([zirc, chrom], [ri, ro, ro+h1], N); % create waveguide description 
guw = cyl.fullyCoupled(n); % waves of circumferential order n
ff = computeW(guw, k)/2/pi; kk = k.*ones(size(ff));
figure, plot(kk(:), ff(:), '.'); ylim([0, 18e3]/(r(end)-r(1)));
xlabel('wavenumber k in rad/m'), ylabel('frequency f in Hz')

% second deposited layer thickness
cyl = Cylinder([zirc, chrom], [ri, ro, ro+h2], N); % create waveguide description 
guw = cyl.fullyCoupled(n); % waves of circumferential order n
ff = computeW(guw, k)/2/pi; kk = k.*ones(size(ff));
hold on, plot(kk(:), ff(:), '.'); ylim([0, 18e3]/(r(end)-r(1)));
legend(strcat(num2str([h1; h2]/1e-6), ' um'), 'Location', 'southeast');