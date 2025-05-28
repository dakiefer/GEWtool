%% Compute leaky Lamb waves in a plate loaded by one fluid
% 
% 
% 2025 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France

mat = MaterialIsotropic('brass'); % load from database (or create your own)
fa = MaterialFluid('water'); 
fe = MaterialIsotropic('pmma');
fb = fa; fb.cl = 1200; fb.rho = 800; 
h = 1e-3;                         % thickness in m
N = 20;                           % number of nodes (dictates accuracy)
w = 2*pi*linspace(1e-3, 7, 500).'*1e6; % frequencies where to compute wavenumbers k
plate = PlateLeaky({mat fe}, [h inf], N);         % create waveguide description 
gew = plate.Lamb; tic;         % choose S+A Lamb waves (assembles matrices)
% gew.op = opExpandTerm(gew, 'Rtop', gew.halfSpaces(2).mat);
gew.op = opExpandTerm(gew, 'Rbottom', gew.halfSpaces(1).mat);
% linearizeInK(gews);             % optional: this makes the computation faster
dat = computeK(gew, w, 120); toc;     % solve 

%%
figure(1); hold on
ph = plot3(real(dat(1).k(:))/1e3, imag(dat(1).k(:))/1e3, dat(1).w(:)/2/pi/1e6, '.','DisplayName','S'); 
% plot3(real(dat(2).k(:))/1e3, imag(dat(2).k(:))/1e3, dat(2).w(:)/2/pi/1e6, '.','DisplayName','A'); 
xlim([0, 12]), ylim([-10.5, 10.5]), view(22, 18)
xlabel('Re(k) in rad/mm'), ylabel('Im(k) in rad/mm'), zlabel('f in MHz')
legend(legendUnq, 'Location', 'southeast')
title(sprintf('Lamb waves in %gmm %s', h/1e-3, mat.name))
