%% Plot the modal field of a bilayered anisotropic waveguide
% The displacements and the y-tractions are visualized. These need to be
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
T = stress(gew, datMode);       % nodal stress

figure(1), % dispersion curves
ph = plot(dat.k(indk,indw), dat.w(indk,indw)/2/pi, 'rd'); % show selected in disp. curve
legend(ph, {'selection'}, 'Location', 'southeast')

% interploate onto finer grid for plotting purposes:
[ui, yi] = GEWinterpolate(gew, u, 200); % interpolate onto 200 equi-distant coordinates
uix = ui(:,1); uiy = ui(:,2);           % interpolated displacement components
Ti = GEWinterpolate(gew, T, yi);        % interpolate stress on same grid
Tiyx = Ti(:,2,1); Tiyy = Ti(:,2,2);     % interpolated stress components

%%  plot the modal displacements 
% All displacement components need to be continuous across the layer interfaces.
figure(2), clf

subplot(2,1,1); hold on, % plot ux
plot(real(uix), yi/1e-3, 'DisplayName', 'Re $u_x$');
plot(imag(uix), yi/1e-3, 'DisplayName', 'Im $u_x$');
yline(gew.geom.yItf(2)/1e-3, '-', {mat2.name, mat1.name}, ...
    'LineWidth', 1,'LabelVerticalAlignment', 'middle', 'HandleVisibility','off')
ylim([gew.geom.yItf(1), gew.geom.yItf(end)]/1e-3)
xlabel('modal displacement ux'), ylabel('y in mm')
legend('Location','northeastoutside')

subplot(2,1,2); hold on % plot uy
plot(real(uiy), yi/1e-3, 'DisplayName', 'Re $u_y$');
plot(imag(uiy), yi/1e-3, 'DisplayName', 'Im $u_y$');
yline(gew.geom.yItf(2)/1e-3, '-', {mat2.name, mat1.name}, ...
    'LineWidth', 1,'LabelVerticalAlignment', 'middle', 'HandleVisibility','off')
ylim([gew.geom.yItf(1), gew.geom.yItf(end)]/1e-3)
xlabel('modal displacement uy'), ylabel('y in mm')
legend('Location','northeastoutside')

%% plot the modal stresses
% The stress components Tyx and Tyy need to be continuous across the layers. Moreover, 
% they vanish at the boundaries. 
figure(3), clf,

subplot(2,1,1); hold on, % plot Tyx
plot(real(Tiyx), yi/1e-3, 'DisplayName', 'Re $T_{yx}$');
plot(imag(Tiyx), yi/1e-3, 'DisplayName', 'Im $T_{yx}$');
yline(gew.geom.yItf(2)/1e-3, '-', {mat2.name, mat1.name}, ...
    'LineWidth', 1,'LabelVerticalAlignment', 'middle', 'HandleVisibility','off')
ylim([gew.geom.yItf(1), gew.geom.yItf(end)]/1e-3)
xlabel('modal stress Tyx'), ylabel('y in mm')
legend('Location','northeastoutside')

subplot(2,1,2); hold on, % plot Tyy
plot(real(Tiyy), yi/1e-3, 'DisplayName', 'Re $T_{yy}$');
plot(imag(Tiyy), yi/1e-3, 'DisplayName', 'Im $T_{yy}$');
yline(gew.geom.yItf(2)/1e-3, '-', {mat2.name, mat1.name}, ...
    'LineWidth', 1,'LabelVerticalAlignment', 'middle', 'HandleVisibility','off')
ylim([gew.geom.yItf(1), gew.geom.yItf(end)]/1e-3)
xlabel('modal stress Tyy'), ylabel('y in mm')
legend('Location','northeastoutside')
