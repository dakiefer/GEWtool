% Test if the universal anisotropy index AU is computed correctly. 
% Also test the Voigt and Reuss homogenization.
% 
% Run using: runtests(). see also: 
% https://fr.mathworks.com/help/matlab/matlab_prog/write-script-based-unit-tests.html
%
% 2024 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France

%% test isotropic material
mat0 = Material('aluminum'); 
assert(mat0.AU == 0);
mat = mat0.rotateEuler(pi/7,0,0); % turn should not affect AU
assert(mat.AU == 0);
mat = mat0.rotateEuler(0,0,pi/7); 
assert(mat.AU == 0);

%% test silicon (cubic)
AU = 0.244; % computed by formula (rounded to 3 digits after the comma)
mat0 = Material('silicon'); 
assert(mat0.AU == AU);
mat = mat0.rotateEuler(pi/7,0,0); % turn to strange off-axis direction (homogenization instead of formula)
assert(mat.AU == AU);

%% test T800/913 (transversely isotropic)
AU = 15.795; % we computed this once and think this is correct  (rounded to 3 digits after the comma)
mat0 = Material('T800_913'); 
assert(mat0.AU == AU);
mat = mat0.rotateEuler(pi/7,0,0); 
assert(mat.AU == AU);
mat = mat0.rotateEuler(0,pi/5,0); 
assert(mat.AU == AU);
mat = mat0.rotateEuler(0,0,pi/3); 
assert(mat.AU == AU);

%% test triclinic
AU = 4.760; % we computed this once and think this is correct  (rounded to 3 digits after the comma)
mat0 = Material('triclinic'); 
assert(mat0.AU == AU);
mat = mat0.rotateEuler(pi/7,0,0); 
assert(mat.AU == AU);
mat = mat0.rotateEuler(0,pi/5,0); 
assert(mat.AU == AU);
mat = mat0.rotateEuler(0,0,pi/3); 
assert(mat.AU == AU);
