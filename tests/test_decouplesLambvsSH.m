% Run using: runtests()
% Test the testing on decoupling of Lamb (flexural/longitudinal) and SH
% (torsional) polarizations.
%
% see also: 
% https://fr.mathworks.com/help/matlab/matlab_prog/write-script-based-unit-tests.html
%
% 2023 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris | PSL, France

r = [10-3.56, 10]*1e-3; 
h = r(end) - r(1);
N = 5;
mat = Material('silicon');
matRot = mat.rotateEuler(22.5/180*pi,'z');

%% test plate 
plate = Plate(mat,h,N);
decoupl = plate.decouplesLambvsSH;
assert(decoupl == true)
plate = Plate(matRot,h,N);
decoupl = plate.decouplesLambvsSH;
assert(decoupl == false);

%% test cylinder 
% zeroth-order circumferential waves
n = 0;
cyl = Cylinder(mat,r,N);
decoupl = cyl.decouplesLambvsSH(n);
assert(decoupl == true);
cyl = Cylinder(matRot,r,N);
decoupl = cyl.decouplesLambvsSH(n);
assert(decoupl == false);

% first-order circumferential waves
n = 1;
cyl = Cylinder(mat,r,N);
decoupl = cyl.decouplesLambvsSH(n);
assert(decoupl == false);
cyl = Cylinder(matRot,r,N);
decoupl = cyl.decouplesLambvsSH(n);
assert(decoupl == false);
