%% Dispersion calculation in a piezoelectric plate
% The material file contains the piezoelectric stress constants "e" (in Voigt
% notation "E") and is loaded with the class "MaterialPiezoelectric". When the
% resulting variable "mat" is passed to the Plate class, GEWtool automatically
% accounts for the piezoelectric coupling. By default, the plate is electrically
% open on both surfaces.
% 
% _↑z_________________________________________ top surface
%  →x   --> k  guided wave    c,e,epsilon,rho  plate
% ____________________________________________ bottom surface
% 
% stiffness c, piezoelectric const. e, permittivity epsilon, mass density rho
%
% See also:
% D. A. Kiefer, G. Watzl, K. Burgholzer, M. Ryzy, and C. Grünsteidl,
% “Electroelastic guided wave dispersion in piezoelectric plates: spectral
% methods and laser-ultrasound experiments,” Nov. 2024. doi:
% 10.48550/arXiv.2412.07389.
% 
% 2024-2025 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France

mat = MaterialPiezoelectric('lithium_niobate'); % Z-cut, X-propagation
h = 1e-3;                        % thickness in m
N = 12;                          % number of nodes (dictates accuracy)
k = linspace(1e-2, 15, 200)/h;   % wavenumbers to solve for
plate = Plate(mat, h, N);        % create waveguide description 
gew = plate.fullyCoupled; tic;   % full polarization

%% plate with electrically open surfaces (no metallization, the default)
datOpen = computeW(gew, k, 20); toc; % solve for 20 modes

figure(1); clf; hold on;
title(sprintf('%s, %g mm thick', mat.name, h/1e-3))
plot(datOpen);

%% plate with electrically shorted surfaces (metallized)
dofPotential = gew.geom.gdofBC{1}(4,:); % the 4th-component are the potenetials, top and bottom
gew = gew.fixGdof(dofPotential); tic    % set potential at boundaries to zero
datShorted = computeW(gew, k, 20); toc;

plot(datShorted, ':')
legend({'elec. open', 'elec. shorted'})
