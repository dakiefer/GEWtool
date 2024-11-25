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

installdir = fileparts(fileparts(which('test_path')));
addpath('.');
addpath(fullfile(installdir, 'GEWdat'));
addpath(fullfile(installdir, 'layers'));
addpath(fullfile(installdir, 'material'));
addpath(fullfile(installdir, 'material/database'));
addpath(fullfile(installdir, 'resourcesAndDeps'));
addpath(fullfile(installdir, 'solvers'));
addpath(fullfile(installdir, 'utilitiesAndNumerics'));
addpath(fullfile(installdir, 'waveguides'));

% common variables:
h = 1e-3;
N = 8;
k = linspace(1e-2, 12, 2)/h; % wavenumbers to solve for

%% Plate example:
mat = Material('steel');
plate = Plate(mat, h, N);
gew = plate.Lamb;
dat = computeW(gew, k); 
cg = energyVelAxial(gew,dat); % check postprocessing

% % If no errors were thrown, we assume that everything is on the path as it
% should.

%% restore path
% This is run even if the previous test fail (when using "runtests")
restoredefaultpath; 
addpath(pp);
