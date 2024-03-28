% Run using: runtests()
% Test cylinder: single layer vs. double layer of same material and thickness.
% Before running, you can set 
% >> show = true; 
% to visualize the output.
%
% see also: 
% https://fr.mathworks.com/help/matlab/matlab_prog/write-script-based-unit-tests.html
%
% 2022 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France

% specify parameters:
r = [20e-3, 20.5e-3, 21e-3]; % radial interface coordinates of layers in m
N = [15, 15]; % number of discretization points
mat = MaterialIsotropic('steel');
n = 0; % circumferential order
nModes = 10;
k = linspace(1e-2, 15, 200)/(r(end)-r(1)); % wavenumber (solve for frequency)
relTol = 1e-6;

%% one vs two layers
% % one layer:
cyl = Cylinder(mat, [r(1), r(3)], sum(N));
gew = cyl.fullyCoupled(n);
dat0 = computeW(gew, k, nModes); 

% % two layers of same material and same total thickness as one layer:
cyl = Cylinder([mat, mat], r, N);
gew = cyl.fullyCoupled(n);
dat = computeW(gew, k, nModes); 

err = abs(dat.w - dat0.w)./dat0.w;
assert( all(err < relTol, "all") );

% % plot 
if exist('show', 'var') && show
    figure, plot(dat0.k(:), dat0.w(:)/2/pi, 'gx'); ylim([0, 4e3]/(r(end)-r(1)));
    xlabel('wavenumber k in rad/m'), ylabel('frequency f in Hz')
    hold on, plot(dat.k(:), dat.w(:)/2/pi, 'k.'); ylim([0, 4e3]/(r(end)-r(1)));
    xlabel('wavenumber k in rad/m'), ylabel('frequency f in Hz')
    legend({'single', 'two lay.'}, 'Location','southeast')
end

