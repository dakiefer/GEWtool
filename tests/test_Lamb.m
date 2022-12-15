% run using: runtests()
% test the Lamb wave solutions vs. precomputed reference values
%
% 2022 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France

% common variables to all tests:
load('data/reference_Lamb.mat') % loads: datRef, h, mat, N
steel = Material('steel');
k = linspace(1e-1, 12, 200)/h; % wavenumbers to solve for
nModes = size(datRef.w, 2); % number of modes
w0 = median(datRef.w, 'all'); % for normalization
assert( mat == steel ) % fails if material has changed

%% single layer - frequency
plate = Plate(steel, h, 40); % uses smaller discretization
gew = plate.Lamb;
dat = computeW(gew, k, nModes); 
Deltaw = abs(dat.w - datRef.w)/w0; % normalized to w0 (median)
assert( all(Deltaw < 1e-8, 'all') )

%% two vs one layer - frequency
hs = [0.71, 0.29]*h; % to test bilayered plate
Ns = [24, 19]; % bilayered plate: number of discretization points

plate = Plate([steel, steel], hs, Ns);
gew = plate.Lamb;
dat = computeW(gew, k, nModes); 
Deltaw = abs(dat.w - datRef.w)/w0; % normalized to w0 (median)
assert( all(Deltaw < 1e-8, 'all') )

% % plot single layer vs. doulbe layer
figure, plot(datRef.k(:), datRef.w(:)/2/pi, 'x'); ylim([0, 6e3]/(sum(hs)));
xlabel('wavenumber k in rad/m'), ylabel('frequency f in Hz')
hold on, plot(dat.k(:), dat.w(:)/2/pi, '.');  ylim([0, 6e3]/(sum(hs)));
legend({'single', 'two lay.'}, 'Location','southeast')
