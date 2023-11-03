%% Dispersion surface of an anisotropic silicon plate close to the ZGV points
% Computes and plots the dispersion surface of the S1/S2b modes close to the
% zero-group-velocity (ZGV) points. The result is plotted as iso-frequency contours.
%
% Literature: 
% D. A. Kiefer, S. Mezil, and C. Prada, “Beating resonance patterns and
% orthogonal wave propagation due to zero-group-velocity guided elastic waves.”
% arXiv, Jul. 26, 2023. doi: 10.48550/arXiv.2307.14259.
% 
% 2023 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France

mat0  = Material('silicon'); % load from database
mat0  = mat0.rotateEuler(0, 45/180*pi, 0);   % x-axis aligned with [110] crystal axis
mat45 = mat0.rotateEuler(0, -45/180*pi, 0);  % x-axis aligned with [010] crystal axis
h = 524.6e-6;   % plate thickness 
N = 16;         % number of nodes (discretization)
k = linspace(2.4, 4.2, 50)*1e3; % wavenumber list to solve for

%% compute ZGV points along both principal axes
% compute ZGV at  0°/[110] (ZGV2, saddle point): 
plate = Plate(mat0, h, N);
gew = plate.fullyCoupledS;    % symmetric waves only (Lamb and SH coupled)
zgv0 = computeZGVScan(gew);   % compute ZGV point (frequency and wavenumber)
w0 = zgv0.w(1); k0 = zgv0.k(1); 
% compute ZGV at 45°/[010] (ZGV1, minimum):
plate = Plate(mat45, h, N);
gew = plate.fullyCoupledS;    % symmetric waves only (Lamb and SH coupled)
zgv45 = computeZGVScan(gew);  % compute ZGV point (frequency and wavenumber)
w45 = zgv45.w(1); k45 = zgv45.k(1);
% combine ZGV data for later plotting:
thetaZGV0 = [0,pi/2,pi,-pi/2];
thetaZGV45 = thetaZGV0+pi/4;
kX_ZGV = [k0*cos(thetaZGV0), k45*cos(thetaZGV45)];
kY_ZGV = [k0*sin(thetaZGV0), k45*sin(thetaZGV45)];

%% compute dispersion surface (in cylindrical coordinates):
thetas = 2*pi*linspace(0, 1, 240); % compute dispersion curves at these angles
clear dat; clear tgv; tic; % clear pre-existing data structures, if any
for i = 1:length(thetas)
    mat = mat0.rotateEuler(0, -thetas(i), 0); % rotate material 
    plate = Plate(mat, h, N);
    gew = plate.fullyCoupledS;        % symmetric waves only (Lamb and SH coupled)
    dat(i) = computeW(gew, k, 3);     % compute the first 3 modes
    tgv(i) = computeZGV(gew, w0, k0); % compute TGV waves with initial guess (w0,k0)
end
toc; % print computational time
kX_TGV = [tgv.k].*cos(thetas); % convert to Cartisian coordinates
kY_TGV = [tgv.k].*sin(thetas); % convert to Cartisian coordinates

% prepare for plotting:
[kk, tt] = meshgrid(k, thetas); % plotting function needs a grid
kX=kk.*cos(tt); % convert to Cartisian coordinates
kY=kk.*sin(tt); % convert to Cartisian coordinates
W=zeros(size(kX)); % allocate for frequencies
nMode = 3; % index of mode to plot
for i = 1:length(thetas), W(i,:)=dat(i).w(:,nMode); end

%% contour plot
figure(34); clf; hold on, grid off; axis equal;
title('S1/S2b mode dispersion contours')
plot([-4 4 nan -4 4 nan -4 4 nan 0 0],[-4 4 nan 4 -4 nan 0 0 nan 4 -4],...
    'LineWidth',1,'Color',[.4,.4,.4],'HandleVisibility','off'); % reference lines
levels = 1/2/pi/1e6*linspace(w45*(1+1e-4), w0*(1 + 1.5e-3), 8); % iso-frequency list
contour(kX/1e3,kY/1e3,W/2/pi/1e6,levels,'LineWidth',2,'DisplayName','iso-freq.'); % contours
plot(kX_TGV/1e3, kY_TGV/1e3, 'r', 'LineWidth', 3,'DisplayName','TGV');     % TGV waves
plot(kX_ZGV/1e3, kY_ZGV/1e3, '.r', 'MarkerSize', 40, 'DisplayName','ZGV'); % ZGV resonances
text(1, 0.35, "[110]"); % crystallographic axis [110]
text(0.5, 1, "[010]", "Rotation",45); % crystallographic axis [010]
cb = colorbar; cb.Label.String = "frequency $\omega/2\pi$ in MHz";
cb.Label.Interpreter = "latex";
legend('Location','southeast');
xlim([-1, 1]*max(kX_TGV)/1e3*1.2), ylim([-1, 1]*max(kY_TGV)/1e3*1.2) % set limits
clim([0.998*w45/2/pi/1e6, 1.001*w0/2/pi/1e6]);
xlabel('$k_{X}$ in rad/mm','Interpreter','latex')
ylabel('$k_{Y}$ in rad/mm','Interpreter','latex')
zlabel('$\omega/2\pi$ in MHz','Interpreter','latex')
