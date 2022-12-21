% Run using: runtests()
% Test if the path contains all needed functions.
%
% see also: 
% https://fr.mathworks.com/help/matlab/matlab_prog/write-script-based-unit-tests.html
%
% 2022 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France

% setup clean path:
pp = path; % save current path
restoredefaultpath; 

installdir = fullfile('../..'); % adjust to path where you put GEWtool!
addpath(fullfile(installdir, 'GEWtool', 'GEWdat'));
addpath(fullfile(installdir, 'GEWtool', 'layers'));
addpath(fullfile(installdir, 'GEWtool', 'material'));
addpath(fullfile(installdir, 'GEWtool', 'material/database'));
addpath(fullfile(installdir, 'GEWtool', 'numerics'));
addpath(fullfile(installdir, 'GEWtool', 'resourcesAndDeps'));
addpath(fullfile(installdir, 'GEWtool', 'solvers'));
addpath(fullfile(installdir, 'GEWtool', 'waveguides'));

% common variables:
h = 1e-3;
N = 12;
k = linspace(1e-2, 12, 10)/h; % wavenumbers to solve for

%% Plate example:
mat = Material('steel');
plate = Plate(mat, h, N);
gew = plate.Lamb;
dat = computeW(gew, k); 

% % If no errors were thrown, we assume that everything is on the path as it
% should.

%% restore path
% This is run even if the previous test fail (when using "runtests")
addpath(pp);
