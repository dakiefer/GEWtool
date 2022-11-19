%% compute guided ultrasonic waves in a layered plate
% 
% 2022 Daniel Kiefer
% Institut Langevin, Paris, France
% 

steel = Material('steel');
chrom = Material('chromium');
hs = [0.71e-3, 0.29e-3]; % thickness list
Ns = [15, 12]; % number of discretization points
k = linspace(1e-2, 12, 200)/(sum(hs)); % wavenumbers to solve for

%% one layer:
plate = Plate(steel, sum(hs), sum(Ns));
gew = plate.Lamb;
dat = computeW(gew, k); 
figure, plot(dat.k(:), dat.w(:)/2/pi, 'x'); ylim([0, 6e3]/(sum(hs)));
xlabel('wavenumber k in rad/m'), ylabel('frequency f in Hz')

%% two layers of same material and same total thickness as one layer:
plate = Plate([steel, steel], hs, Ns);
gew = plate.Lamb;
dat = computeW(gew, k); 
hold on, plot(dat.k(:), dat.w(:)/2/pi, '.');  ylim([0, 6e3]/(sum(hs)));
xlabel('wavenumber k in rad/m'), ylabel('frequency f in Hz')
legend({'single', 'two lay.'}, 'Location','southeast')
