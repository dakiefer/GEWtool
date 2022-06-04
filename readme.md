# GEWTool

Dispersion curves and computation with guided elastic waves (GEWs) in MATLAB©. 

## Installation 

1. Add the folders `layers`, `material`, `material/database`, `numerics`, `solvers` and `waveguides` to the matlab path and save for future sessions.  You could use:
  ```matlab
  installdir = fullfile('~', 'Documents', 'MATLAB'); % adjust to path where you put GEWTool!
  addpath(fullfile(installdir, 'GEWTool', 'layers'));
  addpath(fullfile(installdir, 'GEWTool', 'material'));
  addpath(fullfile(installdir, 'GEWTool', 'material/database'));
  addpath(fullfile(installdir, 'GEWTool', 'numerics'));
  addpath(fullfile(installdir, 'GEWTool', 'solvers'));
  addpath(fullfile(installdir, 'GEWTool', 'waveguides'));
  savepath % make permanent
  ```
2. Enjoy!

## Disclaimer 

2022 Daniel A. Kiefer, [daniel.kiefer@espci.fr](mailto:daniel.kiefer@espci.fr)
Institut Langevin, ESPCI Paris, Paris, France

The functions `collocD`, `fclencurt` and `lglnodes` created by Greg von Winckel are bundled together with their license file and are also available at https://fr.mathworks.com/matlabcentral/profile/authors/869721.