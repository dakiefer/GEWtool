%% Compute leaky Lamb waves in a plate loaded by one fluid
% 
% 
% 2025 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France

mat = MaterialIsotropic('plateMat', 81.3120e9, 40.6560e9, 8.4e3); % brass
matA = MaterialIsotropic('soft', 2.6785e9, 0.6655e9, 2.2e3); % teflon
% mat = MaterialIsotropic('brass'); % load from database (or create your own)
% matA = MaterialFluid('water'); 
% fe = MaterialIsotropic('pmma');
% fb = fa; fb.cl = 1200; fb.rho = 800; 
h = 1e-3;                         % thickness in m
N = 14;                           % number of nodes (dictates accuracy)
% w = 2*pi*(0:25e3:7e6).'; % as in the paper
w = 2*pi*linspace(1e-3, 7, 100).'*1e6; % frequencies where to compute wavenumbers k
plate = PlateLeaky({ matA mat matA}, [ inf h inf], N);         % create waveguide description 
% plate.np.c0 = 1; plate.np.rho0 = 1; plate.np.fh0 = 1; 
gew = plate.Lamb;         % choose S+A Lamb waves (assembles matrices)
opts.subspace = true; 
opts.show = true;
opts.standardEVP = false; % NOT WORKING
tic;
dat = computeK(gew, w, 260, opts); 
time = toc;     % solve 
fprintf("time per frequency is %gs\n",time/length(w));


%% compare to previous
% % plot
addpath('~/Projekte/Radiation/leaky wave calc');
folder = fileparts(which('getMatricesFEMsolid2'));
fileName = fullfile(folder, 'results','oldResultsSBFEM','plate_BrassTeflon','compare_cp_fewer_points.fig');
open(fileName); hold on;
% indRef = abs(imag(kRef)) < 0.5*abs(real(kRef)); % & real(betah) > 1e-4 & real(etah) > 1e-4;
% plot(ffRef(indRef), 2*pi*ffRef(indRef)./real(kRef(indRef)), 'ob', 'MarkerSize',6,'DisplayName','multipar');
ind = imag(dat.k*h) > 1e-3 & imag(dat.k*h) < 0.5*abs(real(dat.k*h)); % & real(betah) < -1e-4 & real(etah) < -1e-4; % is empty!
plot(dat.w(ind)/2/pi/1e6, dat.w(ind)./real(dat.k(ind))/1e3, 'xg', 'MarkerSize',6,'LineWidth',2,'DisplayName','state');
legend; drawnow;
% % plot attenuation
fileName = fullfile(folder, 'results','oldResultsSBFEM','plate_BrassTeflon','compare_att_fewer_points.fig');
open(fileName);
hold all
% att=imag(kRef)/h*20/log(10)*1000;
% plot(ffRef(indRef), att(indRef), 'ob', 'MarkerSize',6,'DisplayName','MultiPar');
att=imag(dat.k)*20/log(10);
plot(dat.w(ind)/2/pi/1e6, att(ind), 'xg', 'MarkerSize',6,'LineWidth',2,'DisplayName','state');
drawnow;

%%
% figure; hold on
% ph = plot3(real(dat(1).k(:))/1e3, imag(dat(1).k(:))/1e3, dat(1).w(:)/2/pi/1e6, '.','DisplayName','S'); 
% % plot3(real(dat(2).k(:))/1e3, imag(dat(2).k(:))/1e3, dat(2).w(:)/2/pi/1e6, '.','DisplayName','A'); 
% xlim([0, 12]), ylim([-10.5, 10.5]), view(22, 18)
% xlabel('Re(k) in rad/mm'), ylabel('Im(k) in rad/mm'), zlabel('f in MHz')
% legend(legendUnq, 'Location', 'southeast')
% title(sprintf('Lamb waves in %gmm %s', h/1e-3, mat.name))