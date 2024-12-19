clear variables
close all

%% Material
mat = Material('LiNbO3');
% 128 deg Y-X cut orientation
mat = mat.rotateEuler(0*pi/180, -(128-90)*pi/180, 90*pi/180);


%% Parameters
cR = 3992; % Expected rayleigh wave speed [m/s] (for high frequency . thickness number)
fmin = 0.1e6; % [Hz]
fmax = 15e6; % [Hz]
h = 0.5*1e-3; % thickness of plate

N = 40; % number of nodes (dictates accuracy)
k = linspace(2*pi*fmin/cR, 2*pi*fmax/cR, 2e3); % wavenumbers to solve for
nModes = 8; % nn of modes saved

plate = Plate(mat, h, N);
gews = plate.fullyCoupled;
dat = computeW(gews, k, nModes);

%% Mode selection (k1,w1)
modeN = 1; % 1st Antisymmetric (A) mode
f1 = 4*1e6; % [Hz]
w1 = 2*pi*f1; % [rad/s]
k1 = get_k(dat.w(:,modeN), dat.k(:,modeN), f1); % [rad/m]
c1 = w1/k1; % [m/s]

% frequency . thickness number is low
% phase speed will be < cR
fprintf('f * h = %0.1f km/s \n',  f1 * h/1e3)
fprintf('cR = %0.2f m/s \n', cR)
fprintf('At f = %0.1f MHz, phase velocity of mode 1 = %0.2f m/s \n', f1/1e6, c1)
fprintf('relative ratio (c1-cR)/cR = %0.2f %% \n', 100*abs(cR-c1)/cR)

%% plot frequency-wavenumber dispersion
figPos = [1921,28,1920,963];
matName = strrep(mat.name,'_',' ');

fig = figure('Position',figPos,'Units','pixels'); clf; hold on
plot(dat.k/1e3, dat.w/2/pi/1e6, '.-');
% xlim([0, 60]);
ylim([0, 2*f1/1e6])
xlabel('wavenumber k in rad/mm'), ylabel('frequency f in MHz')
title(sprintf('Dispersion curve in %gmm %s', h/1e-3, matName))
xline(k1/1e3,'--',sprintf('%0.2f rad/mm',k1/1e3),'FontSize',18)
yline(f1/1e6,'--',sprintf('%d MHz',f1/1e6),'FontSize',18)
set(gca,'Color','none','FontName','Arial','FontSize',18,'FontWeight',...
    'bold','LineWidth',1);

%% plot phase velocity-frequency dispersion
figure('Position',figPos,'Units','normalized'); clf; hold on;
plot(dat.w/2/pi/1e6, dat.w./dat.k/1e3, '.-');
xlim([0, 2*f1/1e6])
ylim([0, 13000/1e3])
xlabel('frequency f in MHz'), ylabel('phase velocity in mm/us')
title(sprintf('Phase velocity curve in %gmm %s', h/1e-3, matName))
xline(f1/1e6,'--',sprintf('%d MHz',f1/1e6),'FontSize',18)
yline(c1/1e3,'--',sprintf('%0.2f mm/us',c1/1e3),'FontSize',18)
yline(cR/1e3,'--',sprintf('cR = %0.2f mm/us',cR/1e3),'FontSize',18)
set(gca,'Color','none','FontName','Arial','FontSize',18,'FontWeight',...
    'bold','LineWidth',1);

%% Check angle dependance on phase velocity of first mode at (w1,k1)
% phase velocity on 128' rot Y-LiNbO3 vs. angle from X-axis

% Angles for rotation along X axis
thetas = 0:1:90;

fmin = 0.1e6; % [Hz]
fmax = 8e6; % [Hz]
N = 20; % number of nodes (dictates accuracy)
k = linspace(2*pi*fmin/cR, 2*pi*fmax/cR, 1e3); % wavenumbers to solve for
nModes = 3; % nn of modes saved

cp = zeros(1, length(thetas));
for i = 1:length(thetas)
    mat = Material('LiNbO3');
    mat = mat.rotateEuler(0*pi/180, -38*pi/180, 90*pi/180);

    % rotate along X axis
    mat = mat.rotateEuler(thetas(i)*pi/180, 'x');

    % Compute waves in a plate
    plate = Plate(mat, h, N);
    gews = plate.fullyCoupled;
    dat = computeW(gews, k, nModes);

    % estimation of phase velocity for mode 1
    kModeN = get_k(dat.w(:,modeN), dat.k(:,modeN), f1); % [rad/m]
    cp(i) = w1/kModeN; % [m/s]
end

figure;
plot(thetas,cp);
xlim([0 90]);
% ylim([0 cR]);
xlabel('angle (deg)'); ylabel('phase velocity [m/s]')
title(sprintf('Mode %d', modeN))

fprintf('phase velocity = %0.2f +- %0.2f m/s \n', mean(cp), std(cp))
fprintf('relative ratio max min cp = %0.2f %% \n', 100 * abs(max(cp) - min(cp)) / max(cp))

%%
function k0 = get_k(w, k, f0)
x = w;
y = k;
x_interp = 2*pi*f0;
k0 = interp1(x, y, x_interp, "linear");
end
