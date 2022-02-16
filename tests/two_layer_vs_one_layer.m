% specify parameters:
r = [4e-3, 4.995e-3, 5e-3];
N = [15, 15];
Nudof = [3, 3];
zirc = Material('zircaloy', 99.3e9, 0.37, 6560, 'Enu'); 
steel = Material('steel');
mats = [zirc, steel];
n = 0;
kh = linspace(1e-3, 5, 700); % wavenumber-thickness (solve for frequency)

%% only first material
c0 = mats(1).c(1,2,1,2); h0 = r(3) - r(1); % normalization parameters
rho0 = mats(1).rho; f0 = sqrt(c0/rho0)/h0; % normalization parameters
geom = Geometry([r(1), r(3)], N(1), Nudof);
l1 = LayerCylindrical(mats(1), r(1), r(3), N(1));
[M, L0, L1, L2] = assembleLayers(geom, l1, n, c0, h0, rho0);
[M, L0, L1, L2] = freeBCs(M, L0, L1, L2, geom, l1, n, c0, h0);

% solve for frequency:
whn = solveWithKh(kh, M, L0, L1, L2);
fh = real(whn/2/pi*f0*h0); 

% plot wavenumbers:
kkh = kh.*ones(size(fh));
hold on, plot(kkh(:)/h0, fh(:)/h0, 'gx');

%% only second material
c0 = mats(2).c(1,2,1,2); h0 = r(3) - r(1); % normalization parameters
rho0 = mats(2).rho; f0 = sqrt(c0/rho0)/h0; % normalization parameters
geom = Geometry([r(1), r(3)], N(2), Nudof);
l1 = LayerCylindrical(mats(2), r(1), r(3), N(2));
[M, L0, L1, L2] = assembleLayers(geom, l1, n, c0, h0, rho0);
[M, L0, L1, L2] = freeBCs(M, L0, L1, L2, geom, l1, n, c0, h0);

% solve for frequency:
whn = solveWithKh(kh, M, L0, L1, L2);
fh = real(whn/2/pi*f0*h0); 

% plot wavenumbers:
kkh = kh.*ones(size(fh));
hold on, plot(kkh(:)/h0, fh(:)/h0, 'cx');

%% bilayer problem thick-thin:
kh = linspace(1e-3, 5, 400); % wavenumber-thickness (solve for frequency)
b = linspace(r(1)*(1+1e-4), r(2)*(1-1e-4), 20);
cc = inferno(length(b));
for ii = 1:length(b)
    r0 = [4e-3, b(ii), 5e-3];
    c0 = mats(1).c(1,2,1,2); h0 = (r0(end)-r0(1))/length(mats); % normalization parameters
    rho0 = mats(1).rho; f0 = sqrt(c0/rho0)/h0; % normalization parameters

    geom = Geometry(r0, N, Nudof);
    l1 = LayerCylindrical(mats(1), r0(1), r0(2), N(1));
    l2 = LayerCylindrical(mats(2), r0(2), r0(3), N(2));
    [M, L0, L1, L2] = assembleLayers(geom, [l1, l2], n, c0, h0, rho0);
    [M, L0, L1, L2] = freeBCs(M, L0, L1, L2, geom, [l1, l2], n, c0, h0);

    % solve for frequency:
    whn = solveWithKh(kh, M, L0, L1, L2);
    fh = real(whn/2/pi*f0*h0); 

    % plot wavenumbers:
    kkh = kh.*ones(size(fh));
    hold on, plot(kkh(:)/h0, fh(:)/h0, '.', 'Color', cc(ii,:));
end 

xlim([0, 5e3]), ylim([0, 2e6])
ylabel('f in Hz'), xlabel('k in rad/m')
legend({'zircaloy', 'steel'})

