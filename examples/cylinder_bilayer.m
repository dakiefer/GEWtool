%% compute axially guided ultrasonic waves in bilayered cylinder
% Waves in a cylinder with two layeres of different material.
% Two different deposited layer thicknesses are compared.
% 
% 2022 Daniel Kiefer
% Institut Langevin, Paris, France
% 

clear all
%% specify parameters:
ri = 5e-3; ro = 5.57e-3; % inner and outer radii
h1 = 10e-6; h2 = 30e-6; % two different deposited layer thicknesses 
N = [20, 10]; % number of discretization points
n = 0; % circumferential order
k = linspace(1e-3, 20, 300)/(ro - ri); % wavenumbers to solve for (use instead of freq)
freq = linspace(1e-3, 14, 300).'*1e6;  % frequencies to solve for (use instead of k)
zirc = Material('zircaloy');
chrom = Material('chromium');

%% with thin layer:
cyl = Cylinder([zirc, chrom], [ri, ro, ro+h1], N); % create waveguide description 
guw = cyl.fullyCoupled(n); % waves of circumferential order n
dat = computeW(guw, k); 
figure, plot(dat.k(:), dat.w(:)/2/pi, '.'); ylim([0, 4e3]/guw.h0);
xlabel('wavenumber k in rad/m'), ylabel('frequency f in Hz')

%% with thick layer:
cyl = Cylinder([zirc, chrom], [ri, ro, ro+h2], N); % create waveguide description 
guw = cyl.fullyCoupled(n); % waves of circumferential order n
dat = computeW(guw, k); 
hold on, plot(dat.k(:), dat.w(:)/2/pi, '.'); 
xlim([0, 20/(ro-ri)]); 
ylim([0, 4e3]/guw.h0); 
legend(strcat(num2str([h1; h2]/1e-6), ' um'), 'Location', 'southeast');
