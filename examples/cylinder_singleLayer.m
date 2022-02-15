%% compute axially guided ultrasonic waves in a cylinder
% different circumferential orders n of the waves are compared
% 
% 2022 Daniel Kiefer
% Institut Langevin, Paris, France
% 

% specify parameters:
steel = Material('steel'); % material
r = [15e-3, 16e-3]; h = r(end) - r(1); % radii and thickness 
N = 15; % number of discretization points
k = linspace(1e-2, 12, 200)/h; % wavenumbers to solve for
cyl = Cylinder(steel, r, N); % create waveguide description 

% compute and plot
n = 0:1:5; % circumferential order to solve for
figure, hold on; ax = gca;
xlabel('wavenumber k in rad/m'), ylabel('frequency f in Hz')
ax.ColorOrder = parula(length(n));
for n0 = n
    guw = cyl.fullyCoupled(n0); % waves of circumferential order n0
    ff = computeW(guw, k)/2/pi; kk = k.*ones(size(ff));
    plot(kk(:), ff(:), '.'); ylim([0, 6e3]/h); drawnow;
end
legend(strcat('n = ', num2str(n.')), 'Location', 'southeast')
