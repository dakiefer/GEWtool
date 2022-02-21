%% compare to zheng
% axial waves in a cylinder of an anisotropic steel material. 
%
% see: M. Zheng, C. He, Y. Lyu, and B. Wu, “Guided waves propagation in anisotropic 
% hollow cylinders by Legendre polynomial solution based on state-vector formalism,” 
% Composite Structures, vol. 207, pp. 645–657, Jan. 2019, 
% doi: 10.1016/j.compstruct.2018.09.042.


b = 141.3/2*1e-3; % outer radius
h = 12.5e-3; % thickness
a = b-h; % inner radius
N = 20; % number of collocation points
% fh = linspace(1e-2, 6, 300).'*1e3;
k = linspace(1e-2, 15, 200)/h;
mat = Material('steel_zheng');
cyl = Cylinder(mat, [a, b], N);

% plot reference:
load('data/zheng.mat');
figure, hold on, ylim([0, 12e3]), xlim([0, 0.4])
plot(L01.f, L01.cp, 'k'), plot(L02.f, L02.cp, 'k'), plot(L04.f, L04.cp, 'k')
plot(F11.f, F11.cp, 'k'), plot(F12.f, F12.cp, 'k'), plot(F13.f, F13.cp, 'k')
xlabel('f in MHz'), ylabel('cp in m/s'), drawnow
ax = gca; ax.ColorOrderIndex = 1; % reset color order 
legend({'', '', '', '', '', 'Zheng'})

% zeroth order:
n = 0; % circumferential wavenumber (flexural order)
guw = cyl.fullyCoupled(n);
ff = computeW(guw, k)/2/pi; kk = k.*ones(size(ff));
plot(ff(:)/1e6, 2*pi*ff(:)./kk(:), 'r.'); 
legend({'', '', '', '', '', 'Zheng', 'n = 0'})

% first order:
n = 1; % circumferential wavenumber (flexural order)
guw = cyl.fullyCoupled(n);
ff = computeW(guw, k)/2/pi; kk = k.*ones(size(ff));
plot(ff(:)/1e6, 2*pi*ff(:)./kk(:), 'g.'); 
legend({'', '', '', '', '', 'Zheng', 'n = 0', 'n = 1'})
