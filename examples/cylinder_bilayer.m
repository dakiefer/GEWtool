%% compute axially guided ultrasonic waves in bilayered cylinder
% Waves in a cylinder with two layeres of different material.
% Two different deposited layer thicknesses are compared.
% 
% 2022 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France

% % specify parameters:
ri = 5e-3; ro = 5.57e-3; % inner and outer radii
h1 = 10e-6; h2 = 30e-6; % two different deposited layer thicknesses 
N = [20, 10]; % number of discretization points
n = 0; % circumferential order
k = linspace(1e-3, 20, 150)/(ro - ri); % wavenumbers to solve for (use instead of freq)
zirc = Material('zircaloy');
chrom = Material('chromium');

%% with thin layer:
cyl = Cylinder([zirc, chrom], [ri, ro, ro+h1], N); % create waveguide description 
gew = cyl.fullyCoupled(n); % waves of circumferential order n
dat = computeW(gew, k); 
cc = lines(2);
figure, ph1 = plot(dat.k, dat.w/2/pi, 'Color', cc(1,:)); ylim([0, 4e3]/gew.np.h0);
xlabel('wavenumber k in rad/m'), ylabel('frequency f in Hz')

%% with thick layer:
cyl = Cylinder([zirc, chrom], [ri, ro, ro+h2], N); % create waveguide description 
gew = cyl.fullyCoupled(n); % waves of circumferential order n
dat = computeW(gew, k); 
hold on, ph2 = plot(dat.k, dat.w/2/pi, 'Color', cc(2,:)); 
xlim([0, 20/(ro-ri)]); 
ylim([0, 4e3]/gew.np.h0); 
legend([ph1(1), ph2(1)], strcat(num2str([h1; h2]/1e-6), ' um'), 'Location', 'southeast');
