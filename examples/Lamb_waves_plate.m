%% Compute Lamb waves in a plate
% This example shows how to compute Lamb waves. The first approach shown
% computes all Lamb waves. The second block computes the symmetric and
% anti-symmetric waves separately. This is faster and the modes are correctly
% orderd. The third section computes the complex wavenumber spectrum. 
% 
% 2022-2023 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France

mat = MaterialIsotropic('steel');        % load from database (or create your own)
h = 1e-3;                       % thickness in m
N = 14;                         % number of nodes (dictates accuracy)
k = linspace(1e-2, 12, 100)/h;  % wavenumbers to solve for
plate = Plate(mat, h, N);       % create waveguide description 

%% compute Lamb waves (any symmetry):
gew = plate.Lamb; tic           % choose waves (assemble matrices) 
dat = computeW(gew, k); toc     % solve
figure(1); clf; 
plot(dat.k/1e3, dat.w/2/pi/1e6); ylim([0, 6]);   % plot
xlabel('wavenumber k in rad/mm'), ylabel('frequency f in MHz')
title('Lamb waves')

%% compute symmetric/anti-symmetric waves separately:
gews = plate.LambSA; tic; % assembles matrices for both the sym/anti-sym Lamb waves
datSA = computeW(gews, k); toc;
figure(2); clf; hold on, cc = lines(3);
hS = plot(datSA(1).k/1e3, datSA(1).w/2/pi/1e6, 'Color', cc(1,:)); % symmetric waves
hA = plot(datSA(2).k/1e3, datSA(2).w/2/pi/1e6, 'Color', cc(2,:)); % anti-symmetric waves
ylim([0, 6]);
xlabel('wavenumber k in rad/mm'), ylabel('frequency f in MHz')
legend([hS(1), hA(1)], {'symmetric', 'anti-symmetric'}, 'Location', 'southeast')
title('symmetric and anti-symmetric Lamb waves')

%% compute wavenumbers and plot:
w = 2*pi*linspace(1e-2, 7, 1000).'*1e6; tic; % frequencies where to compute wavenumbers k
% linearizeInK(gews); % optional: this makes the computation faster
datCompl = computeK(gews, w); toc;
figure(3); clf; hold on, cc = lines(3);
plot3(real(datCompl(1).k(:))/1e3, imag(datCompl(1).k(:))/1e3, datCompl(1).w(:)/2/pi/1e6, '.', 'Color', cc(1,:)); 
plot3(real(datCompl(2).k(:))/1e3, imag(datCompl(2).k(:))/1e3, datCompl(2).w(:)/2/pi/1e6, '.', 'Color', cc(2,:)); 
xlim([0, 12]), ylim([-10.5, 10.5]), view(22, 18)
xlabel('Re(k) in rad/mm'), ylabel('Im(k) in rad/mm')
zlabel('f in MHz')
title('complex spectrum')
