% Test if the universal anisotropy index AU is computed correctly. 
% Also test the Voigt and Reuss homogenization.
% 
% Run using: runtests(). see also: 
% https://fr.mathworks.com/help/matlab/matlab_prog/write-script-based-unit-tests.html
%
% 2024 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France

%% test homogenization
% Verify that homogenization with a uniform orientation distribution leads to an
% isotropic material. 
mat = Material('triclinic');
[voigt,reuss] = homogenizeUniform(mat,10000);
assert(voigt.AU < 1e-3);
assert(reuss.AU < 1e-3);
Cv = voigt.C; Cv = Cv/norm(Cv); Cv = round(Cv,3);
assert(testIsotropy(Cv));
Cr = reuss.C; Cr = Cr/norm(Cr); Cr = round(Cr,3);
assert(testIsotropy(Cr))

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

function ret = testIsotropy(C)
    ret = false; 
    if ~issymmetric(C), return; end
    indZero = [4 5 6 10 11 12 16 17 18 19 20 21]; 
    if C(indZero) ~= 0, return; end
    lbd = C(1,2); 
    if ~(C(1,3) == lbd && C(2,3) == lbd), return; end
    mu = C(4,4); 
    if ~(C(5,5) == mu && C(6,6) == mu), return; end
    if ~(   abs(C(1,1) -lbd-2*mu) < 1e-2 &&...
            abs(C(2,2) -lbd-2*mu) < 1e-2 && ...
            abs(C(3,3) -lbd-2*mu) < 1e-2 ) % use bigger tolerance than in rounding of C
        return; 
    end
    ret = true; 
end