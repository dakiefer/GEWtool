% specify parameters:
clear all
r = [20e-3, 20.5e-3, 21e-3]; % radial interface coordinates of layers in m
N = [15, 15]; % number of discretization points
mat.rho = 7900; lbd = 1.1538e11; mu = 7.6923e10; % steel material
II = eye(3).*shiftdim(eye(3), -2); % 4th order "unit tensor"
mat.c = lbd*II + mu*(permute(II, [1 3 4 2]) + permute(II, [1 3 2 4])); % stiffness tensor
n = 0; % circumferential order
k = linspace(1e-2, 15, 200)/(r(end)-r(1)); % wavenumber (solve for frequency)

%% one vs two layers
% % one layer:
cyl = Cylinder(mat, [r(1), r(3)], sum(N));
gew = cyl.fullyCoupled(n);
dat = computeW(gew, k); 
fig = figure;
plot(dat.k(:), dat.w(:)/2/pi, 'gx'); ylim([0, 4e3]/(r(end)-r(1)));
xlabel('wavenumber k in rad/m'), ylabel('frequency f in Hz')

% % two layers of same material and same total thickness as one layer:
cyl = Cylinder([mat, mat], r, N);
gew = cyl.fullyCoupled(n);
dat = computeW(gew, k); 
hold on, plot(dat.k(:), dat.w(:)/2/pi, 'k.'); ylim([0, 4e3]/(r(end)-r(1)));
xlabel('wavenumber k in rad/m'), ylabel('frequency f in Hz')
legend({'single', 'two lay.'}, 'Location','southeast')

% request user to evaluate test:
assert( userTestConfirmation(fig) )
