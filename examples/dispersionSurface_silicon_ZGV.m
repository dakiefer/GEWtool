%% Dispersion surface of an anisotropic silicon plate close to the ZGV points
% Computes and plots the dispersion surface of the S1/S2b modes close to the
% zero-group-velocity (ZGV) points. The result is plotted as iso-frequency contours.
%
% Literature: 
% D. A. Kiefer, S. Mezil, and C. Prada, "Beating resonance patterns and extreme
% power flux skewing in anisotropic elastic plates," Science Advances, vol. 9,
% no. 51, p. eadk6846, Dec. 2023, doi: 10.1126/sciadv.adk6846.
% 
% 2023 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France

mat0  = Material('silicon'); % load from database
mat0  = mat0.rotateEuler( 45/180*pi, 'z');  % x-axis aligned with [110] crystal axis
mat45 = mat0.rotateEuler(-45/180*pi, 'z');  % x-axis aligned with [010] crystal axis
h = 524.6e-6;   % plate thickness 
N = 16;         % number of nodes (discretization)
nMode = 3;      % index of mode to plot
wmax = 2*pi*10e6;  % maximum angular frequency of interest (for faster computation)

%% compute ZGV points along both principal axes of silicon
% Note that this block will only be valid for cubic materials. adapt if necessary.
% compute ZGV at  0°/[110] (ZGV2, saddle point): 
plate = Plate(mat0, h, N);
gew = plate.fullyCoupledS;        % symmetric waves only (Lamb and SH coupled)
zgv0 = computeZGVScan(gew,wmax);  % compute ZGV point (frequency and wavenumber)
w0 = zgv0.w(1); k0 = zgv0.k(1); 
clear tgv; % clear pre-existing data structures, if any 
tgv(1) = computeZGV(gew, w0, k0); % later we compute all TGVs
% compute ZGV at 45°/[010] (ZGV1, minimum):
plate = Plate(mat45, h, N);
gew = plate.fullyCoupledS;        % symmetric waves only (Lamb and SH coupled)
zgv45 = computeZGVScan(gew,wmax); % compute ZGV point (frequency and wavenumber)
w45 = zgv45.w(1); k45 = zgv45.k(1);
% combine ZGV data for later plotting:
thetaZGV0 = [0,pi/2,pi,-pi/2];
thetaZGV45 = thetaZGV0+pi/4;
kX_ZGV = [k0*cos(thetaZGV0), k45*cos(thetaZGV45)];
kY_ZGV = [k0*sin(thetaZGV0), k45*sin(thetaZGV45)];

%% compute TGV waves in all directions:
thetas = 2*pi*linspace(0, 1, 180); tic; % compute at these angles
for i = 2:length(thetas) % omit first
    mat = mat0.rotateEuler(-thetas(i), 'z');
    plate = Plate(mat, h, N);
    gew = plate.fullyCoupledS; % symmetric waves only (Lamb and SH coupled)
    tgv(i) = computeZGV(gew,tgv(i-1).w,tgv(i-1).k); % compute TGV, initial guess is last solution
end
toc; % print computational time
wZGVmin = min([tgv.w]); wZGVmax = max([tgv.w]);
Dw = wZGVmax-wZGVmin;
kTGVmin = min([tgv.k]); kTGVmax = max([tgv.k]);
kX_TGV = [tgv.k].*cos(thetas); % convert to Cartesian coordinates
kY_TGV = [tgv.k].*sin(thetas); % convert to Cartesian coordinates

%% compute dispersion surface (in cylindrical coordinates) in the TGV region:
k = linspace(kTGVmin*0.1, kTGVmax*1.5, 80); % wavenumber list to solve for
clear dat; tic; % clear pre-existing data structures, if any
for i = 1:length(thetas) % at same angles as before
    mat = mat0.rotateEuler(-thetas(i), 'z');
    plate = Plate(mat, h, N);
    gew = plate.fullyCoupledS; % symmetric waves only (Lamb and SH coupled)
    dat(i) = computeW(gew, k, 3); % compute the first 3 modes
end
toc; % print computational time

% prepare for plotting:
[kk, tt] = meshgrid(k, thetas); % plotting function needs a grid
kX=kk.*cos(tt); % convert to Cartesian coordinates
kY=kk.*sin(tt); % convert to Cartesian coordinates
W=zeros(size(kX)); % allocate for frequencies
for i = 1:length(thetas), W(i,:)=dat(i).w(:,nMode); end

%% contour plot
figure(34); clf; hold on, grid off; axis equal;
title(sprintf('%s (%s)', mat0.name, mat0.symmetry))
klim = max([kX_TGV, kY_TGV])*1.3; 
xlim([-1, 1]*klim/1e3), ylim([-1, 1]*klim/1e3) % set limits
clim([0.998*wZGVmin/2/pi/1e6, 1.001*wZGVmax/2/pi/1e6]);
levels = 1/2/pi/1e6*linspace(wZGVmin*(1+1e-4), wZGVmax+0.3*Dw, 8); % iso-frequency list
contour(kX/1e3,kY/1e3,W/2/pi/1e6,levels,'LineWidth',2,'DisplayName','iso-freq.'); % contours
plot(kX_TGV/1e3, kY_TGV/1e3, 'r', 'LineWidth', 3,'DisplayName','TGV');     % TGV waves
if mat0.symmetry == "cubic"
    plot(klim/1e3*[-1 1 nan -1 1 nan -1 1 nan 0 0],klim/1e3*[-1 1 nan 1 -1 nan 0 0 nan 1 -1],...
        'LineWidth',1,'Color',[.4,.4,.4],'HandleVisibility','off'); % reference lines
    plot(kX_ZGV/1e3, kY_ZGV/1e3, '.r', 'MarkerSize', 40, 'DisplayName','ZGV'); % ZGV resonances
    text(1, 0.35, "[110]"); % crystallographic axis [110]
    text(0.5, 1, "[010]", "Rotation",45); % crystallographic axis [010]
end
cb = colorbar; cb.Label.String = "frequency $\omega/2\pi$ in MHz"; cb.Label.Interpreter = "latex";
legend('Location','southeast','NumColumns',3);
xlabel('$k_{X}$ in rad/mm','Interpreter','latex')
ylabel('$k_{Y}$ in rad/mm','Interpreter','latex')
zlabel('$\omega/2\pi$ in MHz','Interpreter','latex')
