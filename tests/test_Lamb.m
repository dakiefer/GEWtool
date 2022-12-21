% Run using: runtests()
% Test the Lamb wave solutions vs. precomputed reference values.
%
% see also: 
% https://fr.mathworks.com/help/matlab/matlab_prog/write-script-based-unit-tests.html
%
% 2022 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France

% common variables for all tests:
load('data/Lamb_ref.mat') % load pre-solved solutions
steel = Material('steel');
k = kRef; % wavenumbers to solve for
nModes = size(datRef.w, 2); % number of modes
w0 = median(datRef.w, 'all'); % for normalization
assert( matRef == steel ) % fails if material has changed

%% single layer - frequency
plate = Plate(steel, hRef, 40); % uses smaller discretization
gew = plate.Lamb;
dat = computeW(gew, k, nModes); 
Deltaw = abs(dat.w - datRef.w)/w0; % normalized to w0 (median)
assert( all(Deltaw < 1e-8, 'all') )

%% two vs one layer - frequency
hs = [0.71, 0.29]*hRef; % to test bilayered plate
Ns = [24, 19]; % bilayered plate: number of discretization points

plate = Plate([steel, steel], hs, Ns);
gew = plate.Lamb;
dat2 = computeW(gew, k, nModes); 
Deltaw2 = abs(dat2.w - datRef.w)/w0; % normalized to w0 (median)
assert( all(Deltaw2 < 1e-8, 'all') )

% % plot single layer vs. doulbe layer
if exist('show', 'var') && show
    figure, plot(datRef.k(:), datRef.w(:)/2/pi, 'x'); ylim([0, 6e3]/(sum(hs)));
    xlabel('wavenumber k in rad/m'), ylabel('frequency f in Hz')
    hold on, plot(dat2.k(:), dat2.w(:)/2/pi, '.');  ylim([0, 6e3]/(sum(hs)));
    legend({'single', 'two lay.'}, 'Location','southeast')
end
