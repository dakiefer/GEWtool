% % Test multilayer composite plate
% Run using: runtests()
%
% This is a symmetric plate with several layers that is reduced to the top half layers
% for the computation. We have one ZGV point for the symmetric Lamb waves. The
% example is from 
% A. M. A. Huber and M. G. R. Sause, "Classification of solutions
% for guided waves in anisotropic composites with large numbers of layers," The
% Journal of the Acoustical Society of America, vol. 144, no. 6, pp. 3236–3251,
% Dec. 2018, doi: 10.1121/1.5082299.
% 
% 2023 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France

% setup the test case:
mat = Material('T800_913');
p00 =   mat; p00.name = 'p00'; % at 0 degree
p90 =  mat.rotateEuler(0,  90/180*pi, 0); p90.name = 'p90';
p45 =  mat.rotateEuler(0,  45/180*pi, 0); p45.name = 'p45';
m45 = mat.rotateEuler(0, -45/180*pi, 0); m45.name = 'm45';
ply = [p00,p90,p45,m45]; % plys or unit cell that will repeat
mats = [repmat(ply,1,5), fliplr(repmat(ply,1,5))]; % repeat 
hl = 1.25e-3; % layer thickness
Nl = 2*4; % S- and A-waves use half the number of points
plate = Plate(mats, hl*ones(size(mats)), Nl*ones(size(mats)));
h = plate.h;
wmax = 3.5e3*2*pi/h; % maximum frequency of interest (for plotting and ZGV-search)
gew = plate.fullyCoupledS; % Lamb and SH do not decouple

if exist('show', 'var') && show
    ndof = size(gew.op.M,1); fprintf('total number of dofs: %d\n', ndof);
    figure(1); clf; hold on;
    k = linspace(1e-2, 5, 50)/h;
    clear opts; opts.subspace = true; opts.sparse = true; opts.parallel = false; % options for the solver
    dat = computeW(gew, k, 15, opts);
    plot(dat.k/1e3, dat.w/2/pi/1e3, 'k');
    ylim([0, wmax/2/pi/1e3])
    ylabel('frequency $f$ in KHz')
    xlabel('wavenumber $k$ in rad/mm')
    title('layup [0/90/45/-45]5s')
end

%% compute ZGVs
clear opts; % keep in case you execute as a script
if exist('show', 'var') && show, opts.show=true;  else, opts.show=false; end 
opts.kEnd = 2.5/h; % restrict search domain for faster computations
tic, zgv = computeZGVScan(gew, wmax, opts); timing = toc;
nZGV = length(zgv.w(zgv.w < wmax));

if exist('show', 'var') && show
    fprintf('Computed %d ZGV points in %g s.\n', nZGV, timing); 
    figure(1), hold on;
    plot(zgv.k/1e3, zgv.w/2/pi/1e3, 'r*');
    if isfield(zgv, 'k0s'), xline(zgv.k0s/1e3); end
end
assert(nZGV == 1, 'Missed ZGV points.');
assert(timing < 8, 'Slow calculation (previously 6 s).');
