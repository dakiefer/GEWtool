%% Plot the modal field of a bilayered anisotropic waveguide
% The displacements and the y-tractions are visualized. These need to be
% continuous across the layer boundaries. 
% 
% 2022 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France

mat1 = Material('zircaloy'); % transversely isotropic
mat2 = Material('chromium'); % isotropic
h1 = 0.7e-3; h2 = 0.3e-3; % thicknesses
N1 = 20; N2 = 10; % number of nodes
k = linspace(1e-2, 10, 100)/(h1+h2); % wavenumbers for dispersion curves

plate = Plate([mat1, mat2], [h1, h2], [N1, N2]); % plate model
gew = plate.Lamb; % assemble matrices for Lamb-polarized waves
dat = computeW(gew,k,10); % compute solutions and keep 10 modes

% plot dispersion curves:
figure(1), clf, hold on 
plot(dat.k, dat.w/2/pi); ylim([0, 6]*1e6)
xlabel('wavenumber k in rad/m'), ylabel('frequency f in Hz')

%% Retrieve modal field:
% select a mode:
indw = 5; % modes are ordered from low to high frequencies at constant k
[~, indk] = min(abs(k - 3000)); % at this specific wavenumber 
figure(1), 
ph = plot(dat.k(indk,indw), dat.w(indk,indw)/2/pi, 'rd'); % show in disp. curve
legend(ph, {'selection'}, 'Location', 'southeast')

% nodal points:
y1 = gew.geom.y{1}; y2 = gew.geom.y{2};  % nodal points of the layers
y = [y1; y2]; 

% retrieve displacements: 
u1 = squeeze(dat.u{1}(indk,indw,:,:)); % displacements of layer 1
u2 = squeeze(dat.u{2}(indk,indw,:,:)); % displacements of layer 2
ux = [u1(:,1); u2(:,1)];
uy = [u1(:,2); u2(:,2)];

% retrieve the stresses (y-traction is ey*T = [Tyx, Tyy]):
T = stress(gew, dat); % provided by GEWdat
T1 = squeeze(T{1}(indk,indw,:,2,:));
T2 = squeeze(T{2}(indk,indw,:,2,:));
Tyx = [T1(:,1); T2(:,1)];
Tyy = [T1(:,2); T2(:,2)];

%%  plot the modal displacements 
% All displacement components need to be continuous across the layer interfaces.
figure(2), clf
subplot(2,1,1); hold on, % plot ux
plot(real(ux), y/1e-3, '-*');
plot(imag(ux), y/1e-3, '-*');
yline(gew.geom.yItf(2)/1e-3, '-', {mat2.name, mat1.name}, 'LineWidth', 1,...
    'LabelVerticalAlignment', 'middle')
ylim([gew.geom.yItf(1), gew.geom.yItf(end)]/1e-3)
xlabel('modal displacement ux'), ylabel('y in mm')
legend({'real ux', 'imag ux'}, 'Location','best')

subplot(2,1,2); hold on % plot uy
plot(real(uy), y/1e-3, '-*');
plot(imag(uy), y/1e-3, '-*');
yline(gew.geom.yItf(2)/1e-3, '-', {mat2.name, mat1.name}, 'LineWidth', 1,...
    'LabelVerticalAlignment', 'middle')
ylim([gew.geom.yItf(1), gew.geom.yItf(end)]/1e-3)
xlabel('modal displacement uy'), ylabel('y in mm')
legend({'real uy', 'imag uy'}, 'Location','best')

%% plot also the modal stresses
% The stress components Tyx and Tyy need to be continuous across the layers. Moreover, 
% they vanish at the boundaries. 
figure(3), clf,
subplot(2,1,1); hold on, % plot Tyx
plot(real(Tyx), y/1e-3, '-*');
plot(imag(Tyx), y/1e-3, '-*');
yline(gew.geom.yItf(2)/1e-3, '-', {mat2.name, mat1.name}, 'LineWidth', 1,...
    'LabelVerticalAlignment', 'middle')
ylim([gew.geom.yItf(1), gew.geom.yItf(end)]/1e-3)
xlabel('modal stress Tyx'), ylabel('y in mm')
legend({'real Tyx', 'imag Tyx'}, 'Location','best')

subplot(2,1,2); hold on, % plot Tyy
plot(real(Tyy), y/1e-3, '-*');
plot(imag(Tyy), y/1e-3, '-*');
yline(gew.geom.yItf(2)/1e-3, '-', {mat2.name, mat1.name}, 'LineWidth', 1,...
    'LabelVerticalAlignment', 'middle')
ylim([gew.geom.yItf(1), gew.geom.yItf(end)]/1e-3)
xlabel('modal stress Tyy'), ylabel('y in mm')
legend({'real Tyy', 'imag Tyy'}, 'Location','best')
