%% compute axially guided ultrasonic waves in a cylinder
% different circumferential orders n of the waves are compared
% 
% 2022 Daniel Kiefer
% Institut Langevin, Paris, France
% 

clear all
%% specify parameters:
rho = 7900; lbd = 1.1538e11; mu = 7.6923e10; % steel material
r = [10e-3, 11e-3]; h = r(end) - r(1); % radii and thickness 
N = 15; % number of discretization points
k = linspace(1e-2, 12, 200)/h; % wavenumbers to solve for

%% material and waveguide description:
II = eye(3).*shiftdim(eye(3), -2); % 4th order "unit tensor"
c = lbd*II + mu*(permute(II, [1 3 4 2]) + permute(II, [1 3 2 4])); % stiffness tensor
steel.c = c; steel.rho = rho; % material 
cyl = Cylinder(steel, r, N); % create waveguide description 

%% compute and plot
n = 0:1:5; % circumferential order to solve for
figure, hold on; ax = gca; xlim([0, k(end)])
xlabel('wavenumber k in rad/m'), ylabel('frequency f in Hz')
ax.ColorOrder = winter(length(n));
for n0 = n
    guw = cyl.fullyCoupled(n0); % assembles matrices for the specified waves
    ff = computeW(guw, k)/2/pi; kk = k.*ones(size(ff));
    plot(kk(:), ff(:), '.'); ylim([0, 6e3]/h); drawnow;
end
legend(strcat('n = ', num2str(n.')), 'Location', 'southeast')
