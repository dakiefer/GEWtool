%% compute axially guided ultrasonic waves in a cylinder
% waves in a cylinder with two layeres of different material
% 
% 2022 Daniel Kiefer
% Institut Langevin, Paris, France
% 

% specify parameters:
r = [5e-3, 5.37e-3 6e-3]; % radial coordinates of layer interfaces
N = [10, 12]; % number of discretization points
n = 0; % circumferential order
k = linspace(1e-2, 15, 400)/(r(end)-r(1)); % wavenumber-thickness (solve for frequency)

% material description:
II = eye(3).*shiftdim(eye(3), -2); % 4th order "unit tensor"
mat = jsondecode(fileread('../material/database/zircaloy_aniso.json'));
zirc.rho = mat.rho;
zirc.c = voigt2tensor(mat.C); % 4th order tensor c created from Voigt notated C
% zirc.c = mat.lambda*II + mat.mu*(permute(II, [1 3 4 2]) + permute(II, [1 3 2 4])); % stiffness tensor
steel.rho = 7900; lbd = 1.1538e11; mu = 7.6923e10; % steel material
steel.c = lbd*II + mu*(permute(II, [1 3 4 2]) + permute(II, [1 3 2 4])); % stiffness tensor

% compute
cyl = Cylinder([steel, zirc], r, N); % create waveguide description 
guw = cyl.fullyCoupled(n); % waves of circumferential order n
ff = computeW(guw, k)/2/pi; kk = k.*ones(size(ff));

% plot
hold on, plot(kk(:), ff(:), '.'); ylim([0, 6e3]/(r(end)-r(1)));
xlabel('wavenumber k in rad/m'), ylabel('frequency f in Hz')
