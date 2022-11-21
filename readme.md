# GEWtool

Dispersion curves and computation with guided elastic waves (GEWs) in MATLAB©. 

## Installation 

1. Add the folders `layers`, `material`, `material/database`, `numerics`, `solvers` and `waveguides` to the matlab path and save for future sessions.  You could use:
  ```matlab
  installdir = fullfile('~', 'Documents', 'MATLAB'); % adjust to path where you put GEWtool!
  addpath(fullfile(installdir, 'GEWtool', 'GEWdat'));
  addpath(fullfile(installdir, 'GEWtool', 'layers'));
  addpath(fullfile(installdir, 'GEWtool', 'material'));
  addpath(fullfile(installdir, 'GEWtool', 'material/database'));
  addpath(fullfile(installdir, 'GEWtool', 'numerics'));
  addpath(fullfile(installdir, 'GEWtool', 'solvers'));
  addpath(fullfile(installdir, 'GEWtool', 'waveguides'));
  savepath % make permanent
  ```
2. Enjoy!

## Disclaimer 

2022 Daniel A. Kiefer, [daniel.kiefer@espci.fr](mailto:daniel.kiefer@espci.fr)
Institut Langevin, ESPCI Paris, Paris, France

The functions `collocD`, `fclencurt` and `lglnodes` created by Greg von Winckel are bundled together with their license file and are also available at https://fr.mathworks.com/matlabcentral/profile/authors/869721.