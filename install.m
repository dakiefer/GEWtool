% Install GEWtool
%
% Execute this script to add GEWtool permanently to your matlab path. You need
% to navigate to the "GEWtool" folder before doing so. This needs to be done
% only once before using GEWtool. You might need to re-execute the script after
% a major update of GEWtool.
% In order to uninstall GEWtool, in the Matlab GUI go to Home -> Set Path. Then
% remove all paths containing "GEWtool" and save. 
% 
% 2022-2023 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France

installdir = pwd;
addpath('.');
addpath(fullfile(installdir, 'GEWdat'));
addpath(fullfile(installdir, 'layers'));
addpath(fullfile(installdir, 'material'));
addpath(fullfile(installdir, 'material/database'));
addpath(fullfile(installdir, 'resourcesAndDeps'));
addpath(fullfile(installdir, 'solvers'));
addpath(fullfile(installdir, 'utilitiesAndNumerics'));
addpath(fullfile(installdir, 'waveguides'));
savepath % make permanent

% OPTIONAL: 
% It is a good idea to collect material files in a folder outside the GEWtool
% folder. In this way, when replacing GEWtool with a new version, you keep your
% material files. For this end, create a folder where to collect your material
% files and then add it to the matlab path: 
% 
% 1. create folder anywhere you like
% 2. add it to the matlab path: 

% folderPath = '~/Documents/MATLAB/materialDB'; % change this to the full path of your folder
% addpath(folderPath);
% savepath; % make permanent
