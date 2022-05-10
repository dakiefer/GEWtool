% specify parameters:
r = [4e-3, 4.995e-3, 5e-3];
N = [12, 12];
Nudof = [3, 3];
zirc = Material('zircaloy'); 
steel = Material('steel');
mats = [zirc, steel];
n = 0;
k = linspace(1e-3, 5, 200)/(r(end)-r(1)); % wavenumber-thickness (solve for frequency)

%% only first material
cyl = Cylinder(mats(1), [r(1), r(3)], N(1));
guw = cyl.fullyCoupled(n);
dat = computeW(guw, k); 
figure, plot(dat.k(:), dat.w(:)/2/pi, 'gx'); 
xlim([0, 5e3]), ylim([0, 2e6]) 
xlabel('wavenumber k in rad/m'), ylabel('frequency f in Hz')

%% only second material
cyl = Cylinder(mats(2), [r(1), r(3)], N(2));
guw = cyl.fullyCoupled(n);
dat = computeW(guw, k); 
hold on, plot(dat.k(:), dat.w(:)/2/pi, 'cx'); drawnow;

%% bilayer problem thick-thin:
b = linspace(r(1)*(1+1e-4), r(3)*(1-1e-4), 20);
cc = inferno(length(b));
for ii = 1:length(b)
    r0 = [r(1), b(ii), r(3)];
    cyl = Cylinder(mats, r0, N);
    guw = cyl.fullyCoupled(n);
    dat = computeW(guw, k); 
    plot(dat.k(:), dat.w(:)/2/pi, '.', 'Color', cc(ii,:)); drawnow;
end 
legend({'zircaloy', 'steel'}, 'Location','southeast')

