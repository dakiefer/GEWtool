% Run using: runtests()
% Compare axial waves in anisotropic cylinder to solutions by Zheng. 
%
% see: M. Zheng, C. He, Y. Lyu, and B. Wu, “Guided waves propagation in anisotropic 
% hollow cylinders by Legendre polynomial solution based on state-vector formalism,” 
% Composite Structures, vol. 207, pp. 645–657, Jan. 2019, 
% doi: 10.1016/j.compstruct.2018.09.042.
%
% see also: 
% https://fr.mathworks.com/help/matlab/matlab_prog/write-script-based-unit-tests.html
%
% 2022 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France
%

b = 141.3/2*1e-3; % outer radius
h = 12.5e-3; % thickness
a = b-h; % inner radius
N = 20; % number of collocation points
% fh = linspace(1e-2, 6, 300).'*1e3;
k = linspace(1e-2, 15, 200)/h;
mat = Material('steel_zheng');
cyl = Cylinder(mat, [a, b], N);

%% plot 
% plot reference:
load('data/zheng.mat');
fig = figure; hold on, ylim([0, 12e3]), xlim([0, 0.4])
plot(L01.f, L01.cp, 'k'), plot(L02.f, L02.cp, 'k'), plot(L04.f, L04.cp, 'k')
plot(F11.f, F11.cp, 'k'), plot(F12.f, F12.cp, 'k'), plot(F13.f, F13.cp, 'k')
xlabel('f in MHz'), ylabel('cp in m/s'), drawnow
ax = gca; ax.ColorOrderIndex = 1; % reset color order 
legend({'', '', '', '', '', 'Zheng'})

% zeroth order:
n = 0; % circumferential wavenumber (flexural order)
gew = cyl.fullyCoupled(n);
dat = computeW(gew, k);
plot(dat.w(:)/2/pi/1e6, dat.w(:)./dat.k(:), 'r.'); 
legend({'', '', '', '', '', 'Zheng', 'n = 0'})

% first order:
n = 1; % circumferential wavenumber (flexural order)
gew = cyl.fullyCoupled(n);
dat = computeW(gew, k);
plot(dat.w(:)/2/pi/1e6, dat.w(:)./dat.k(:), 'g.'); 
legend({'', '', '', '', '', 'Zheng', 'n = 0', 'n = 1'})

% request user to evaluate test:
assert( userTestConfirmation(fig) )
