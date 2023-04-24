% % Test austenitic steel plate
% Run using: runtests()
%
% This is a plate made of orthotropic austenitic steel. It exhibits lots of ZGV 
% points.
% 
% 2023 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France

mat = Material('steel_austenitic'); % orthotropic 
mat = mat.rotateEuler(90/180*pi, 90/180*pi, 0); % this direction is fun
h = 1e-3; % thickness 
N = 18; % number of discretization points
wmax = 2*pi*15e6; % maximum frequency for plotting and ZGV-search
plate = Plate(mat, h, N); % create waveguide description 

%% symmetric Lamb waves
gew = plate.LambS; % assembles matrices for the specified waves
clear opts; % keep in case you execute as a script
if exist('show', 'var') && show, opts.show=true;  else, opts.show=false; end 
opts.kEnd = 15e3;
tic, zgv = computeZGVScan(gew, wmax, opts); timing = toc;
nZGV = length(zgv.w(zgv.w < wmax));  % print number of ZGV points 
if exist('show', 'var') && show, fprintf('Computed %d ZGV points in %g s.\n', nZGV, timing); end 
assert(nZGV == 6, 'Missed ZGV points.');
% assert(timing < 2, 'Slow calculation (previously less than 1s).');

%% anti-symmetric Lamb waves
gew = plate.LambA; % assembles matrices for the specified waves
clear opts; % keep in case you execute as a script
if exist('show', 'var') && show, opts.show=true;  else, opts.show=false; end 
opts.kEnd = 15e3;
tic, zgv = computeZGVScan(gew, wmax, opts); timing = toc;
nZGV = length(zgv.w(zgv.w < wmax));  % print number of ZGV points 
if exist('show', 'var') && show, fprintf('Computed %d ZGV points in %g s.\n', nZGV, timing); end 
assert(nZGV == 6, 'Missed ZGV points.');
% assert(timing < 2, 'Slow calculation (previously less than 1s).');

%% symmetric Lamb waves
% Test if the search is robust to the location of the target wavenumbers. The
% target wavenumbers change when the increment Dk is changed.
% NOTE: change Dk_default if the implementation in computeZGVScan() is changed. 
gew = plate.LambS; % assembles matrices for the specified waves
clear opts; % keep in case you execute as a script
if exist('show', 'var') && show, opts.show=true;  else, opts.show=false; end 
opts.kEnd = 15e3;
DkList = 450:50:900;
for Dk = DkList
    opts.Dk = Dk;
    tic, zgv = computeZGVScan(gew, wmax, opts); timing = toc;
    nZGV = length(zgv.w(zgv.w < wmax));  % print number of ZGV points 
    if exist('show', 'var') && show, fprintf('Computed %d ZGV points in %g s.\n', nZGV, timing); end 
    assert(nZGV == 6, 'Missed ZGV points.');
%     assert(timing < 2, 'Slow calculation (previously less than 1s).');
end

%% Lamb waves
gew = plate.Lamb; % assembles matrices for the specified waves
clear opts; % keep in case you execute as a script
if exist('show', 'var') && show, opts.show=true;  else, opts.show=false; end 
opts.kEnd = 15e3;
tic, zgv = computeZGVScan(gew, wmax, opts); timing = toc;
nZGV = length(zgv.w(zgv.w < wmax));  % print number of ZGV points 

if exist('show', 'var') && show
    fprintf('Computed %d ZGV points in %g s.\n', nZGV, timing);
    figure(1); clf; hold on;
    k = linspace(1e-2, 20, 100)/h;
    dat = computeW(gew, k, 15, opts);
    plot(dat.k/1e3, dat.w/2/pi/1e6, 'k');
    plot(zgv.k/1e3, zgv.w/2/pi/1e6, 'r*');
    if isfield(zgv, 'k0s'), xline(zgv.k0s/1e3); end
    ylim([0, wmax/2/pi/1e6])
    ylabel('frequency $f$ in MHz')
    xlabel('wavenumber $k$ in rad/mm')
    title('orthotropic steel plate')
end
assert(nZGV == 14, 'Missed ZGV points.');
% assert(timing < 8, 'Slow calculation (previously less than 4s).');
