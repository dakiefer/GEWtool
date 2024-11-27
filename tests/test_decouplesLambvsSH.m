% Run using: runtests()
% Test the decoupling of Lamb (flexural/longitudinal) and SH
% (torsional) polarizations.
%
% see also: 
% https://fr.mathworks.com/help/matlab/matlab_prog/write-script-based-unit-tests.html
%
% 2023-2024 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris | PSL, France

r = [10-3.56, 10]*1e-3; 
h = r(end) - r(1);
N = 5;
matIso = Material('steel');
matTri = Material('triclinic');

%% plate on symmetry axis
% test return value of decouplesLambvsSH()
plate = Plate(matIso,h,N);
decoupl = plate.decouplesLambvsSH;
assert(decoupl == true)

% check wheter a warning about decoupling is thrown 
lastwarn('', ''); % reset warnings
gew = plate.Lamb; 
[warnMsg, warnId] = lastwarn(); 
assert(isempty(warnId) && isempty(warnMsg))

%% plate off-symmetry axis
% test return value of decouplesLambvsSH()
plate = Plate(matTri,h,N);
decoupl = plate.decouplesLambvsSH;
assert(decoupl == false);

% check wheter a warning about decoupling is thrown 
lastwarn('', ''); % reset warnings
gew = plate.Lamb; 
[warnMsg, warnId] = lastwarn(); 
assert(strcmp(warnId, 'GEWTOOL:Waveguide:donotdecouple'))
disp('Warning issued on purpose. Testing.')

%% cylinder longitudinal waves on-symmetry
n = 0;

% test return value of decouplesLambvsSH()
cyl = Cylinder(matIso,r,N);
decoupl = cyl.decouplesLambvsSH(n);
assert(decoupl == true);

% check wheter a warning about decoupling is thrown 
lastwarn('', ''); % reset warnings
gew = cyl.longitudinal; 
[warnMsg, warnId] = lastwarn(); 
assert(isempty(warnId) && isempty(warnMsg))

%% cylinder longitudinal waves off-symmetry
n = 0;

% test return value of decouplesLambvsSH()
cyl = Cylinder(matTri,r,N);
decoupl = cyl.decouplesLambvsSH(n);
assert(decoupl == false);

% check wheter a warning about decoupling is thrown 
lastwarn('', ''); % reset warnings
gew = cyl.longitudinal; 
[warnMsg, warnId] = lastwarn(); 
assert(strcmp(warnId, 'GEWTOOL:Waveguide:donotdecouple'))
disp('Warning issued on purpose. Testing.')

%% cylinder flexural waves on-symmetry
n = 1;

% test return value of decouplesLambvsSH()
cyl = Cylinder(matIso,r,N);
decoupl = cyl.decouplesLambvsSH(n);
assert(decoupl == false);

% check wheter a warning about decoupling is thrown 
lastwarn('', ''); % reset warnings
gew = cyl.Lamb(n); % longitudinal() automatically sets n = 0
[warnMsg, warnId] = lastwarn(); 
assert(strcmp(warnId, 'GEWTOOL:Waveguide:donotdecouple'))
disp('Warning issued on purpose. Testing.')

%% cylinder flexural waves off-symmetry
n = 1; 

% test return value of decouplesLambvsSH()
cyl = Cylinder(matTri,r,N);
decoupl = cyl.decouplesLambvsSH(n);
assert(decoupl == false);

% check wheter a warning about decoupling is thrown 
lastwarn('', ''); % reset warnings
gew = cyl.Lamb(n); % longitudinal() automatically sets n = 0
[warnMsg, warnId] = lastwarn(); 
assert(strcmp(warnId, 'GEWTOOL:Waveguide:donotdecouple'))
disp('Warning issued on purpose. Testing.')

%% cylinder circumferential waves on-symmetry

% test return value of decouplesLambvsSH()
cyl = CylinderCircumferential(matIso,r,N);
decoupl = cyl.decouplesLambvsSH;
assert(decoupl == true);

% check wheter a warning about decoupling is thrown 
lastwarn('', ''); % reset warnings
gew = cyl.Lamb; 
[warnMsg, warnId] = lastwarn(); 
assert(isempty(warnId) && isempty(warnMsg))

%% cylinder circumferential waves off-symmetry

% test return value of decouplesLambvsSH()
cyl = CylinderCircumferential(matTri,r,N);
decoupl = cyl.decouplesLambvsSH;
assert(decoupl == false);

% check wheter a warning about decoupling is thrown 
lastwarn('', ''); % reset warnings
gew = cyl.Lamb; 
[warnMsg, warnId] = lastwarn(); 
assert(strcmp(warnId, 'GEWTOOL:Waveguide:donotdecouple'))
disp('Warning issued on purpose. Testing.')

%% call polarization() directly 
udof = [1 3];

plate = Plate(matIso,h,N);
lastwarn('', ''); % reset warnings
gew = plate.polarization(udof,0); 
[warnMsg, warnId] = lastwarn(); 
assert(isempty(warnId) && isempty(warnMsg))

plate = Plate(matTri,h,N);
lastwarn('', ''); % reset warnings
gew = plate.polarization(udof,0); 
[warnMsg, warnId] = lastwarn(); 
assert(strcmp(warnId, 'GEWTOOL:Waveguide:donotdecouple'))
disp('Warning issued on purpose. Testing.')
