% specify parameters:
r = [4e-3, 4.995e-3, 5e-3];
N = [12, 12];
Nudof = [3, 3];
zirc = Material('zircaloy', 99.3e9, 0.37, 6560, 'Enu'); 
steel = Material('steel');
mats = [zirc, steel];
n = 0;
k = linspace(1e-3, 5, 200)/(r(end)-r(1)); % wavenumber-thickness (solve for frequency)

%% only first material
cyl = Cylinder(mats(1), [r(1), r(3)], N(1));
guw = cyl.fullyCoupled(n);
ff = computeW(guw, k)/2/pi; kk = k.*ones(size(ff));
figure, plot(kk(:), ff(:), 'gx'); ylim([0, 4e3]/(r(end)-r(1)));
xlabel('wavenumber k in rad/m'), ylabel('frequency f in Hz')

%% only second material
cyl = Cylinder(mats(2), [r(1), r(3)], N(2));
guw = cyl.fullyCoupled(n);
ff = computeW(guw, k)/2/pi; kk = k.*ones(size(ff));
hold on, plot(kk(:), ff(:), 'cx'); ylim([0, 4e3]/(r(end)-r(1)));
xlabel('wavenumber k in rad/m'), ylabel('frequency f in Hz')

%% bilayer problem thick-thin:
b = linspace(r(1)*(1+1e-4), r(3)*(1-1e-4), 20);
cc = inferno(length(b));
for ii = 1:length(b)
    r0 = [r(1), b(ii), r(3)];
    cyl = Cylinder(mats, r0, N);
    guw = cyl.fullyCoupled(n);
    ff = computeW(guw, k)/2/pi; kk = k.*ones(size(ff));
    hold on, plot(kk(:), ff(:), '.', 'Color', cc(ii,:)); 
    xlabel('wavenumber k in rad/m'), ylabel('frequency f in Hz');
end 

xlim([0, 5e3]), ylim([0, 2e6])
ylabel('frequency f in Hz'), xlabel('wavenumber k in rad/m')
legend({'zircaloy', 'steel'})

