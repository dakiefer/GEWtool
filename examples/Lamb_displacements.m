%% Plot the modal field of a bilayered anisotropic waveguide
% The displacements and the z-tractions are visualized. These need to be
% continuous across the layer boundaries. 
% 
% 2022-2024 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France

% setup bi-layered waveguide:
mat1 = Material('zircaloy');   % anisotropic
mat2 = Material('chromium');    % isotropic
h1 = 0.8e-3; h2 = 0.2e-3;       % thicknesses
N1 = 8; N2 = 6;                % number of nodes
k = linspace(1e-2, 10, 100)/(h1+h2);  % wavenumbers to compute dispersion curves
plate = Plate({mat1, mat2}, [h1, h2], [N1, N2]); % plate model
gew = plate.Lamb;               % assemble matrices for Lamb-polarized waves
dat = computeW(gew,k,8);       % compute solutions 

% plot dispersion curves:
figure(1), clf, hold on 
plot(dat.k, dat.w/2/pi, 'k'); ylim([0, 4]*1e6)
xlabel('wavenumber k in rad/m'), ylabel('frequency f in Hz')

%% Retrieve modal field:
indw = 5; % modes are ordered from low to high frequencies at constant k
[~, indk] = min(abs(k - 3000)); % the closest to this specific wavenumber 
datMode = extractModes(dat,indk,indw); % returns a data structure consistent to ''dat'' but ony with this mode
u = datMode.u;                  % nodal displacements
T = stress(datMode);       % nodal stress

figure(1), % dispersion curves
ph = plot(dat.k(indk,indw), dat.w(indk,indw)/2/pi, 'rd'); % show selected in disp. curve
legend(ph, {'selection'}, 'Location', 'southeast')

% interploate onto finer grid for plotting purposes:
[ui, zi] = GEWinterpolate(gew, u, 200); % interpolate onto 200 equi-distant coordinates
uix = ui(:,1); uiz = ui(:,2);           % interpolated displacement components
Ti = GEWinterpolate(gew, T, zi);        % interpolate stress on same grid
Tizx = Ti(:,2,1); Tizz = Ti(:,2,2);     % interpolated stress components

%%  plot the modal displacements 
% All displacement components need to be continuous across the layer interfaces.
figure(2), clf

subplot(2,1,1); hold on, % plot ux
plot(real(uix), zi/1e-3, 'DisplayName', 'Re $u_x$');
plot(imag(uix), zi/1e-3, 'DisplayName', 'Im $u_x$');
yline(gew.geom.zItf(2)/1e-3, '-', {mat2.name, mat1.name}, ...
    'LineWidth', 1,'LabelVerticalAlignment', 'middle', 'HandleVisibility','off')
ylim([gew.geom.zItf(1), gew.geom.zItf(end)]/1e-3)
xlabel('modal displacement ux'), ylabel('z in mm')
legend('Location','northeastoutside')

subplot(2,1,2); hold on % plot uy
plot(real(uiz), zi/1e-3, 'DisplayName', 'Re $u_z$');
plot(imag(uiz), zi/1e-3, 'DisplayName', 'Im $u_z$');
yline(gew.geom.zItf(2)/1e-3, '-', {mat2.name, mat1.name}, ...
    'LineWidth', 1,'LabelVerticalAlignment', 'middle', 'HandleVisibility','off')
ylim([gew.geom.zItf(1), gew.geom.zItf(end)]/1e-3)
xlabel('modal displacement uz'), ylabel('z in mm')
legend('Location','northeastoutside')

%% plot the modal stresses
% The stress components Tyx and Tyy need to be continuous across the layers. Moreover, 
% they vanish at the boundaries. 
figure(3), clf,

subplot(2,1,1); hold on, % plot Tyx
plot(real(Tizx), zi/1e-3, 'DisplayName', 'Re $T_{zx}$');
plot(imag(Tizx), zi/1e-3, 'DisplayName', 'Im $T_{zx}$');
yline(gew.geom.zItf(2)/1e-3, '-', {mat2.name, mat1.name}, ...
    'LineWidth', 1,'LabelVerticalAlignment', 'middle', 'HandleVisibility','off')
ylim([gew.geom.zItf(1), gew.geom.zItf(end)]/1e-3)
xlabel('modal stress Tzx'), ylabel('z in mm')
legend('Location','northeastoutside')

subplot(2,1,2); hold on, % plot Tyy
plot(real(Tizz), zi/1e-3, 'DisplayName', 'Re $T_{zz}$');
plot(imag(Tizz), zi/1e-3, 'DisplayName', 'Im $T_{zz}$');
yline(gew.geom.zItf(2)/1e-3, '-', {mat2.name, mat1.name}, ...
    'LineWidth', 1,'LabelVerticalAlignment', 'middle', 'HandleVisibility','off')
ylim([gew.geom.zItf(1), gew.geom.zItf(end)]/1e-3)
xlabel('modal stress Tzz'), ylabel('z in mm')
legend('Location','northeastoutside')
