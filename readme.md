# GUWTool

A MATLAB© Toolbox for guided ultrasonic waves (GUWs). 

## Installation 

1. Install `DMSUITE` by Weideman and Reddy:
   [https://mathworks.com/matlabcentral/fileexchange/29-dmsuite](https://mathworks.com/matlabcentral/fileexchange/29-dmsuite)
   
2. Add the folders `layers`, `material`, `material/database`, `solvers` and `waveguides` to the matlab path and save for future sessions.  For example:

  ```matlab
  installdir = fullfile('~', 'Documents', 'MATLAB'); % adjust to path where you put GUWTool!
  addpath(fullfile(installdir, 'GUWTool', 'layers'));
  addpath(fullfile(installdir, 'GUWTool', 'material'));
  addpath(fullfile(installdir, 'GUWTool', 'material/database'));
  addpath(fullfile(installdir, 'GUWTool', 'solvers'));
  addpath(fullfile(installdir, 'GUWTool', 'waveguides'));
  savepath % make permanent
  ```

## Disclaimer 

2022 Daniel A. Kiefer, [daniel.kiefer@espci.fr](mailto:daniel.kiefer@espci.fr)
Institut Langevin, ESPCI Paris, Paris, France
