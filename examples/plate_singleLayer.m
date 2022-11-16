%% compute guided ultrasonic waves in a plate
% 
% 2022 Daniel Kiefer
% Institut Langevin, Paris, France
% 

% % specify parameters:
rho = 7900; lbd = 1.1538e11; mu = 7.6923e10; % steel material
h = 1e-3; % thickness 
N = 15; % number of discretization points
k = linspace(1e-2, 12, 200)/h; % wavenumbers to solve for

% % material and waveguide description:
II = eye(3).*shiftdim(eye(3), -2); % 4th order "unit tensor"
c = lbd*II + mu*(permute(II, [1 3 4 2]) + permute(II, [1 3 2 4])); % stiffness tensor
steel.c = c; steel.rho = rho; % material 
plate = Plate(steel, h, N); % create waveguide description 

%% compute frequencies and plot:
guw = plate.Lamb; % assembles matrices for the specified waves
dat = computeW(guw, k); 
figure, h0 = plot(dat.k(:), dat.w(:)/2/pi, 'ko'); ylim([0, 6e3]/h);
ax = gca; ax.ColorOrderIndex = 1; % reset color order index
xlabel('wavenumber k in rad/m'), ylabel('frequency f in Hz')

%% compute symmetric/anti-symmetric waves separately:
guws = plate.LambSA; % assembles matrices for both the sym/anti-sym Lamb waves
datSA = computeW(guws, k);

datS = datSA(1); % symmetric waves
hold on, h1 = plot(datS.k(:), datS.w(:)/2/pi, '.'); 
datA = datSA(2); % anti-symmetric waves
hold on, h2 = plot(datA.k(:), datA.w(:)/2/pi, '.'); ylim([0, 6e3]/h);
xlabel('wavenumber k in rad/m'), ylabel('frequency f in Hz')
legend([h0, h1, h2], {'coupled (reference)', 'symmetric', 'anti-symmetric'}, 'Location', 'southeast')

%% compute wavenumbers and plot:
freq = linspace(1e-3, 6, 100).'*1e6;
dat = computeK(guw, 2*pi*freq);
hold on, h3 = plot(real(dat.k(:)), dat.w(:)/2/pi, '.'); xlim([0, max(k)])
legend([h0, h1, h2, h3], {'coupled (reference)', 'symmetric', 'anti-symmetric'...
    'complex k'}, 'Location', 'southeast')


