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
k = 10e3;     % choose a wavenumber in rad/m
dat = computeW(gew, k);

% extract the mode to be animated
nW = 3;              % choose a mode index
nK = 1;              % index of the only k that we computed 
w = dat.w(nK,nW);    % angular frequency
Ni = 10*nW + 10;     % number of interpolation points should increase with mode number
[ui,ri] = GEWinterpolate(gew, dat.u{1}(nK,nW,:,:), Ni);
Nphi = 40;           % number of points on circumference
[ue, phi] = GEWextrapolateCircumference(ui, gew.n, Nphi); % extrapolate along angle phi

%% animate modal displacements
% Re-scale displacements and extract components. 
% Note: longitudinal and torsional waves do not have all displacement components.
scale = max(ri)/10; % how strongly the grid is deformed: change if necessary
ue = ue/max(abs(ue(:)))*scale; 
indx = find(gew.displComp == "ux");     % index to extract axial displacements
indphi = find(gew.displComp == "uphi"); % index to extract angular displacements
indr = find(gew.displComp == "ur");     % index to extract radial displacements
if ~isempty(indx),   ux = ue(:,:,indx);     else ux = zeros(size(ue(:,:,1)));   end
if ~isempty(indphi), uphi = ue(:,:,indphi); else uphi = zeros(size(ue(:,:,1))); end
if ~isempty(indr),   ur = ue(:,:,indr);     else ur = zeros(size(ue(:,:,1)));   end

% initial plot
figure('Position', [100, 600, 600, 350]); 
mPlot = mesh(NaN(2), NaN(2), NaN(2), 'EdgeColor','k'); % initialize
title(sprintf('rod: mode %d', nW));
axis manual;  % do not re-scale
axis equal;   % equal scale
axis off;     % do not show
view(-38,30); % camera view

% iterate over time
T = 2*pi/w;     % period
NT = 1;         % number of periods
t = linspace(0, NT*T, NT*80); % time samples
for i = 2:length(t)
    rDeformed   = ri  + real(ur*exp(1i*w*t(i))); 
    phiDeformed = phi + real(uphi./(ri+eps)*exp(1i*w*t(i))); 
    mPlot.XData = rDeformed.*cos(-phiDeformed)/1e-3; % x-phi-r -> X-Y-Z (clock-wise angles)
    mPlot.YData = rDeformed.*sin(-phiDeformed)/1e-3; % x-phi-r -> X-Y-Z (clock-wise angles)
    mPlot.ZData = real(ux*exp(1i*w*t(i)))/1e-3;
    zlim([-1 1]*1.05*scale/1e-3); xlim([-1 1]*1.15); ylim([-1 1]*1.15);
    drawnow; pause(1/25);
end

%% animate the displacement structure
% ux0 = ux(:,1); uphi0 = uphi(:,1); ur0 = ur(:,1); % at phi = 0
% uRange = max(max(abs([ux0, uphi0, ur0]))); % maximum in displacement
% figure; hold on; xlim([-1.05 1.05]*uRange); ylim([r(1), r(end)]/1e-3)
% title(sprintf('rod: mode %d', nW));
% ylabel('radial coor. r in mm'), xlabel('displacement in 1')
% T = 2*pi/w;     % period
% NT = 1;         % number of periods
% t = linspace(0, NT*T, NT*80); % time samples
% for i = 1:length(t)
%     cla;
%     plot(real(ux0*exp(1i*w*t(i))), ri/1e-3, 'DisplayName', 'ux');
%     plot(real(uphi0*exp(1i*w*t(i))), ri/1e-3, 'DisplayName', 'uphi');
%     plot(real(ur0*exp(1i*w*t(i))), ri/1e-3, 'DisplayName', 'ur');
%     legend; drawnow; 
%     pause(1/25);
% end
