%% compute guided waves in a plate
% 
% 2022 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France

mat = Material('steel');        % load from database (or create your own)
h = 1e-3;                       % thickness in m
N = 12;                         % number of nodes (dictates accuracy)
k = linspace(1e-2, 12, 100)/h;  % wavenumbers to solve for
plate = Plate(mat, h, N);       % create waveguide description 

%% compute Lamb waves (any symmetry):
gew = plate.Lamb; tic           % choose waves (assemble matrices) 
dat = computeW(gew, k); toc     % solve
plot(dat.k/1e3, dat.w/2/pi/1e6); ylim([0, 6]);   % plot
xlabel('wavenumber k in rad/mm'), ylabel('frequency f in MHz')
title('Lamb waves')

%% compute symmetric/anti-symmetric waves separately:
gews = plate.LambSA; % assembles matrices for both the sym/anti-sym Lamb waves
datSA = computeW(gews, k);

datS = datSA(1); % symmetric waves
datA = datSA(2); % anti-symmetric waves
figure, hold on, cc = lines(3);
hS = plot(datS.k/1e3, datS.w/2/pi/1e6, 'Color', cc(1,:)); 
hA = plot(datA.k/1e3, datA.w/2/pi/1e6, 'Color', cc(2,:)); ylim([0, 6]);
xlabel('wavenumber k in rad/mm'), ylabel('frequency f in MHz')
legend([hS(1), hA(1)], {'symmetric', 'anti-symmetric'}, 'Location', 'southeast')
title('symmetric and anti-symmetric Lamb waves')

%% compute wavenumbers and plot:
freq = linspace(1e-2, 6, 500).'*1e6; % frequencies where to compute wavenumbers k
dat = computeK(gew, 2*pi*freq);
figure, plot3(real(dat.k(:))/1e3, imag(dat.k(:))/1e3, dat.w(:)/2/pi/1e6, '.'); 
xlim([0, 12]), ylim([-10.5, 10.5]), view(-22, 18)
xlabel('Re(k) in rad/mm'), ylabel('Im(k) in rad/mm')
zlabel('f in MHz')
title('complex spectrum')
