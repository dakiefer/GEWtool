%% Composite plate with many layers
% The layup of the transversely isotropic fiber composite is [0/90/45/-45]50s, i.e., 
% the indicated rotations of the fibers are repeated 50 times and then also repeated 
% symmetrically. This results in 400 layers in total. Each of the layers has a 
% thickness of 0.125 mm.
%
% The example is from:
% [1] A. M. A. Huber and M. G. R. Sause, "Classification of solutions
% for guided waves in anisotropic composites with large numbers of layers," The
% Journal of the Acoustical Society of America, vol. 144, no. 6, pp. 3236–3251,
% Dec. 2018, doi: 10.1121/1.5082299.
% 
% GEWtool computes the dispersion curves in well under one second on a current
% Notebook. Compared to the computational time reported in [1] (101 min.), this
% corresponds to a speedup of roughly factor 10000 (with some uncertainty due to
% different hardware). This speedup comes largely due to GEWtool knowing that
% waves of a same family do not cross, hence, we can take arbitrarily large
% steps in the wavenumber k. We also exploit subspace eigenvalue solvers, sparse
% matrix representation and parallel computing (on 8 cores). The semi-analytical
% approach of GEWtool allows for trivial parallelization as each of the
% eigenvalue problems is independent.
% 
% 2023 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France

% % plot reference extracted from Huber [1]:
figure(1); clf; hold on;
load("../tests/data/dispersionHuber.mat", "dispersionHuber");
plot(dispersionHuber.f, dispersionHuber.cp, 'o', 'MarkerEdgeColor', [.5,.5,.5], 'MarkerSize', 5)
xlim([0, 3.5]), ylim([0, 8]),
xlabel('frequency-thickness $f h$ in MHz mm')
ylabel('phase velocity $c_\mathrm{p}$ in mm/$\mu$s')
title('layup [0/90/45/-45]50s'); 

% % material layup:
mat = Material('T800_913');    % load matrial data from file "T800_913.json" 
% mat = mat.rotateEuler(0, 0/180*pi, 0);
p00 =   mat; p00.name = 'p00'; % fibers at 0 degrees (propagation direction) 
p90 =  mat.rotateEuler(0,  -90/180*pi, 0); p90.name = 'p90'; % fibers at 90 degrees
p45 =  mat.rotateEuler(0,  -45/180*pi, 0); p45.name = 'p45';
m45 = mat.rotateEuler(0, +45/180*pi, 0); m45.name = 'm45';
plys = [p00,p90,p45,m45]; % plys or unit cell that will repeat
mats = [repmat(plys,1,50), fliplr(repmat(plys,1,50))]; % repeat as desired

% % geometry and other parameters:
hl = 0.125e-3; % layer thickness
Nl = 2*(2);    % number of nodes: S- and A-waves use Nl/2 each, minimum is 4/2 = 2.
plate = Plate(mats, hl, Nl); % create geometry and FE matrices.
h = plate.h;   % total thickness
whmax = 3.5e3*2*pi; % maximum frequency of interest (for plotting and ZGV-search)

% % calculate dispersion curves
gew = plate.fullyCoupledSA; % Lamb and SH waves are coupled. Decompose into symmetric/antisymmetric.
fprintf('Total number of dofs for S/A Lamb waves: %d, for %d layers.\n', size(gew(1).op.M,1), plate.geom.nLay);
k = linspace(1e-2, 20, 400)/h; % calculate on these wavenumbers
nModes = 6;           % number of modes to compute
clear opts;           % clear options for the solver (in case opts exists already)
opts.subspace = true; % use eigs() instead of eig(): solve for eigensolutions one by one.
opts.sparse = true;   % use sparse matrices, i.e., zeros are not represented explicitly.
opts.parallel = true; if opts.parallel, gcp(); end % use multicore computing and start a parallel pool
tic, dat = computeW(gew, k, nModes, opts); toc % compute and time (first computation takes a bit longer)

% % plot on top of reference
figure(1), hold on;
plot(dat(1).w*h/2/pi/1e3, dat(1).w./dat(1).k/1e3, 'r'); % symmetric waves
plot(dat(2).w*h/2/pi/1e3, dat(2).w./dat(2).k/1e3, 'b'); % anti-symmetric waves
xlim([0, 3.5]), ylim([0, 8]), drawnow;

%% compute ZGVs
% Our efficient ZGV solver can handle the big matrices that are involved in this
% example. There is only one ZGV in the symmetric waves: 
sel = 1;            % only S-waves exhibit a ZGV point
clear opts,         % clear opts in case it exists already
opts.kEnd = 2.5/h;  % search only up to k*h = 2.5 -> speedup calculation
opts.Neigs = 2;     % significant speed-up compared to default of 5.
opts.show = true;   % show some information on the progress
tic, zgv = computeZGVScan(gew(sel), whmax/h, opts); toc

%%
graphics('paper');
fh = figure(2); clf; hold on;
fh.Position = [50, 50, 600, 450];
ph1 = plot(dat(sel).k*h, dat(sel).w*h/2/pi/1e3, 'k');
ph2 = plot(zgv.k*h, zgv.w*h/2/pi/1e3, 'rd', 'MarkerSize',8,'MarkerFaceColor','r');
ylim([0, whmax/2/pi/1e3])
xlim([0, 5])
legend([ph1(1), ph2], {'S-waves', 'ZGV points'}, 'Location','southeast');
ylabel('frequency-thickness $f h$ in MHz mm')
xlabel('wavenumber-thickness $k h$ in rad')
% title('layup [0/90/45/-45]50s')

