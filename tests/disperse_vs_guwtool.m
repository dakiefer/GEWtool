% parameters:
a = 20; b = 21; h = b-a;
rho=7932; mu=rho*3260^2; lbd=rho*5960^2-2*mu;	
mat = Material('noname', lbd, mu, rho);
N = 12;
k = linspace(1e-3, 0.5, 100)/h; % wavenumber-thickness (solve for frequency)
cyl = Cylinder(mat, [a, b], N);

%% zeroth-order circumferential waves:
n = 0;
guw = cyl.fullyCoupled(n);
dat = computeW(guw, k); 

% plot against reference data:
load('data/disperse.mat')
figure, hold on
plot(disperse.k/1e-3, disperse.f*1e6, 'g.')
xlim([0, 0.5]), ylim([0, 200])
xlabel('wavenumber k in rad/m'), ylabel('frequency f in Hz')
ax = gca; ax.ColorOrderIndex = 1; % reset color order 
hold on, plot(dat.k(:), dat.w(:)/2/pi, '.'); 
legend({'disperse', 'n = 0'})

n = 1;
guw = cyl.fullyCoupled(n);
dat = computeW(guw, k); 
hold on, plot(dat.k(:), dat.w(:)/2/pi, '.'); 
legend({'disperse', 'n = 0'})

%% fist-order circumferential waves:
% n = 1;
% guw = cyl.fullyCoupled(n);
% dat = computeW(guw, k); 
% hold on, plot(dat.k(:), dat.w(:)/2/pi, '.'); 
% legend({'disperse', 'n = 0', 'n = 1'})
