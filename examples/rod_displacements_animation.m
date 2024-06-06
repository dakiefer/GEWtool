%% Animate the modal field of a rod
% The displacements on an arbitrary cross-section are visualized for one selected mode.
% 
% 2024 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France

mat = Material('aluminum');
N = 20;         % number of nodes on radius
h = 1e-3;       % outer radius in m
r = [0, h];     % radial coordinates (full cylinder when inner radius is zero)
guide = Cylinder(mat,r,N); 
gew = guide.longitudinal; 
dat = computeW(gew, 10e3);  % compute at k = 10e3 rad/m

% extract the mode to be animated
nW = 3;              % choose a mode
nK = 1;              % we computed only at one k
w = dat.w(nK,nW);    % angular frequency
Ni = 10*nW + 10;     % number of interpolation points should increase with mode number
[ui,ri] = GEWinterpolate(gew, dat.u{1}(nK,nW,:,:), Ni);
n = 0;               % circumferential order needs to match the wave modeled by ''gew''
Nphi = 40;           % number of points on circumference
[ue, phi] = GEWextrapolateCircumference(ui, n, Nphi); % extrapolate along angle phi

%% animate modal displacements
% re-scale displacements:
scale = max(ri)/20;
ue = ue/max(abs(ue(:)))*scale; 

% initial data at t = 0
rDeformed = ri + real(ue(:,:,2)); 
X = rDeformed.*cos(phi)/1e-3;
Y = rDeformed.*sin(-phi)/1e-3;
Z = real(ue(:,:,1))/1e-3;

% initial plot
figure('Position', [100, 600, 600, 350]); clf; 
mPlot = mesh(X, Y, Z, 'EdgeColor','k');
title(sprintf('rod: mode %d', nW));
zlim([-1.05 1.05]*scale/1e-3);
axis off;

% iterate over time
T = 2*pi/w;     % period
NT = 3;         % number of periods
t = linspace(0, NT*T, NT*80); % time samples
for i = 2:length(t)
    rDeformed = ri + real(ue(:,:,2)*exp(1i*w*t(i))); 
    mPlot.XData = rDeformed.*cos(phi)/1e-3;
    mPlot.YData = rDeformed.*sin(-phi)/1e-3;
    mPlot.ZData = real(ue(:,:,1)*exp(1i*w*t(i)))/1e-3;
    pause(1/25);
end

%% animate the displacement structure
% ux = ui(:,1);
% ur = ui(:,2);
% uRange = max(max(abs([ux, ur]))); % maximum in displacement
% figure; hold on; xlim([-1.05 1.05]*uRange); ylim([r(1), r(end)]/1e-3)
% title(sprintf('rod: mode %d', nW));
% ylabel('radial coor. r in mm'), xlabel('displacement in 1')
% T = 2*pi/w;     % period
% NT = 2;         % number of periods
% t = linspace(0, NT*T, NT*80); % time samples
% for i = 1:length(t)
%     cla;
%     plot(real(ux*exp(1i*w*t(i))), ri/1e-3, 'DisplayName', 'ux');
%     plot(real(ur*exp(1i*w*t(i))), ri/1e-3, 'DisplayName', 'ur');
%     legend; drawnow; 
%     pause(1/25);
% end
