% % Test isotropic plate
% Run using: runtests()
%
% Search the ZGV points in a common isotropic material. 
% 
% 2023 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France

mat = MaterialIsotropic('aluminum');  
h = 1e-3; % thickness 
N = 25; % number of discretization points
wmax = 2*pi*20e6; % maximum frequency for plotting and ZGV-search
plate = Plate(mat, h, N); % create waveguide description 

%% symmetric Lamb waves
gew = plate.LambS; % assembles matrices for the specified waves
clear opts; % keep in case you execute as a script
if exist('show', 'var') && show, opts.show=true;  else, opts.show=false; end 
tic, zgv = computeZGVScan(gew, wmax, opts); timing = toc;
nZGV = length(zgv.w(zgv.w < wmax));  % print number of ZGV points 

if exist('show', 'var') && show
    ndof = size(gew.op.M,1); fprintf('total number of dofs: %d\n', ndof);
    k = linspace(1e-2, 20, 100)/h;
    dat = computeW(gew, k, 15);
    fprintf('Computed %d ZGV points in %g s.\n', nZGV, timing); 
    figure; hold on;
    plot(dat.k/1e3, dat.w/2/pi/1e6, 'k');
    ylim([0, wmax/2/pi/1e6])
    ylabel('frequency $f$ in MHz')
    xlabel('wavenumber $k$ in rad/mm')
    title('symmetric Lamb waves')
    plot(zgv.k/1e3, zgv.w/2/pi/1e6, 'r*');
    if isfield(zgv, 'k0s'), xline(zgv.k0s/1e3); end
end
assert(nZGV == 3, 'Missed ZGV points.');
assert(timing < 4, 'Slow calculation (previously 2 s).');

%% anti-symmetric Lamb waves
gew = plate.LambA; % assembles matrices for the specified waves
clear opts; % keep in case you execute as a script
if exist('show', 'var') && show, opts.show=true;  else, opts.show=false; end 
tic, zgv = computeZGVScan(gew, wmax, opts); timing = toc;
nZGV = length(zgv.w(zgv.w < wmax));  % print number of ZGV points 

if exist('show', 'var') && show
    ndof = size(gew.op.M,1); fprintf('total number of dofs: %d\n', ndof);
    k = linspace(1e-2, 20, 100)/h;
    dat = computeW(gew, k, 15);
    fprintf('Computed %d ZGV points in %g s.\n', nZGV, timing); 
    figure; hold on;
    plot(dat.k/1e3, dat.w/2/pi/1e6, 'k');
    ylim([0, wmax/2/pi/1e6])
    ylabel('frequency $f$ in MHz')
    xlabel('wavenumber $k$ in rad/mm')
    title('anti-symmetric Lamb waves')
    plot(zgv.k/1e3, zgv.w/2/pi/1e6, 'r*');
    if isfield(zgv, 'k0s'), xline(zgv.k0s/1e3); end
end
assert(nZGV == 0, 'Missed ZGV points.');
assert(timing < 10, 'Slow calculation (previously 6 s).');

%% Lamb waves
gew = plate.Lamb; % assembles matrices for the specified waves
clear opts; % keep in case you execute as a script
if exist('show', 'var') && show, opts.show=true;  else, opts.show=false; end 
tic, zgv = computeZGVScan(gew, wmax, opts); timing = toc;
nZGV = length(zgv.w(zgv.w < wmax));  % print number of ZGV points 

if exist('show', 'var') && show
    ndof = size(gew.op.M,1); fprintf('total number of dofs: %d\n', ndof);
    k = linspace(1e-2, 20, 100)/h;
    dat = computeW(gew, k, 20);
    fprintf('Computed %d ZGV points in %g s.\n', nZGV, timing); 
    figure; hold on;
    plot(dat.k/1e3, dat.w/2/pi/1e6, 'k');
    ylim([0, wmax/2/pi/1e6])
    ylabel('frequency $f$ in MHz')
    xlabel('wavenumber $k$ in rad/mm')
    title('Lamb waves')
    plot(zgv.k/1e3, zgv.w/2/pi/1e6, 'r*');
    if isfield(zgv, 'k0s'), xline(zgv.k0s/1e3); end
end
assert(nZGV == 3, 'Missed ZGV points.');
assert(timing < 16, 'Slow calculation (previously 9 s).');
