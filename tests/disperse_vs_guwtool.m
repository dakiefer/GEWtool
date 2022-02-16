% reference data:
load('../reference implementations/cylindrical/comparison_disperse/disperse.mat')
figure, hold on
plot(disperse.k/1e-3, disperse.f*1e6, 'g.')
xlim([0, 0.5]), ylim([0, 200])
xlabel('wavenumber k in rad/m'), ylabel('frequency f in Hz')
ax = gca; ax.ColorOrderIndex = 1; % reset color order 

% parameters:
a = 20; b = 21; h = b-a;
rho=7932; mu=rho*3260^2; lbd=rho*5960^2-2*mu;	
mat = Material('noname', lbd, mu, rho);
N = 25;
k = linspace(1e-3, 0.5, 100)/h; % wavenumber-thickness (solve for frequency)
cyl = Cylinder(mat, [a, b], N);

%% zeroth-order circumferential waves:
n = 0;
guw = cyl.fullyCoupled(n);
ff = computeW(guw, k)/2/pi; kk = k.*ones(size(ff));
hold on, plot(kk(:), ff(:), '.'); 
legend({'disperse', 'n = 0'})

%% fist-order circumferential waves:
n = 1;
guw = cyl.fullyCoupled(n);
ff = computeW(guw, k)/2/pi; kk = k.*ones(size(ff));
hold on, plot(kk(:), ff(:), '.'); 
legend({'disperse', 'n = 0', 'n = 1'})
