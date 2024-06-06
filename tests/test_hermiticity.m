% Test Hermiticity
% Nondissipative waveguides should yield a Hermitian problem. When the problem is not 
% set up carefully, this might be broken due to contamination with numerical noise. 
%
% Run using: runtests(). see also: 
% https://fr.mathworks.com/help/matlab/matlab_prog/write-script-based-unit-tests.html
%
% 2024 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France

matIso = MaterialIsotropic('steel');
matIso = matIso.rotateEuler(pi/5, pi/10, 0); % test if rotation keeps symmetries of tensor
matTri = Material('triclinic'); 
matTri = matTri.rotateEuler(pi/7,pi/10,0);       % use orientation where waves do not decouple
h = 1e-3;
N = 10;

%% test single-layered isotropic plate
guide = Plate(matIso,h,N);
gew = guide.fullyCoupled;
assert(testHermiticity(gew.op));
gew = guide.fullyCoupledS;
assert(testHermiticity(gew.op));
gew = guide.fullyCoupledA;
assert(testHermiticity(gew.op));
gew = guide.Lamb;
assert(testHermiticity(gew.op));
gew = guide.LambS;
assert(testHermiticity(gew.op));
gew = guide.LambA;
assert(testHermiticity(gew.op));
gew = guide.sh;
assert(testHermiticity(gew.op));

%% test single-layered anisotropic plate
guide = Plate(matTri,h,N);
gew = guide.fullyCoupled;
assert(testHermiticity(gew.op));

%% test double-layered isotropic plate
guide = Plate([matIso, matIso], h/2, round(N/2));
gew = guide.fullyCoupled;
assert(testHermiticity(gew.op));
gew = guide.fullyCoupledS;
assert(testHermiticity(gew.op));
gew = guide.fullyCoupledA;
assert(testHermiticity(gew.op));
gew = guide.Lamb;
assert(testHermiticity(gew.op));
gew = guide.LambS;
assert(testHermiticity(gew.op));
gew = guide.LambA;
assert(testHermiticity(gew.op));
gew = guide.sh;
assert(testHermiticity(gew.op));

%% test double-layered anisotropic plate
guide = Plate([matTri, matIso],h,N);
gew = guide.fullyCoupled;
assert(testHermiticity(gew.op));

%% test single-layered isotropic cylinder
guide = Cylinder(matIso,[h, 2*h],N); 
gew = guide.flexural(0);   % n = 0
assert(testHermiticity(gew.op));
gew = guide.longitudinal;  % n = 0
assert(testHermiticity(gew.op));
gew = guide.torsional;     % n = 0
assert(testHermiticity(gew.op));
gew = guide.flexural(1);   % n = 1, generates complex matrices! 
assert(testHermiticity(gew.op));

%% test single-layered anisotropic cylinder
guide = Cylinder(matTri,[h, 2*h],N); 
gew = guide.flexural(0);   % n = 0
assert(testHermiticity(gew.op));
gew = guide.flexural(1);   % n = 1, generates complex matrices! 
assert(testHermiticity(gew.op));

%% test double-layered isotropic cylinder
guide = Cylinder([matIso,matIso],[h, 2*h, 3*h],N); 
gew = guide.flexural(0);   % n = 0
assert(testHermiticity(gew.op));
gew = guide.longitudinal;  % n = 0
assert(testHermiticity(gew.op));
gew = guide.torsional;     % n = 0
assert(testHermiticity(gew.op));
gew = guide.flexural(1);   % n = 1, generates complex matrices! 
assert(testHermiticity(gew.op));

%% test single-layered anisotropic cylinder
guide = Cylinder([matTri,matIso],[h, 2*h, 3*h],N); 
gew = guide.flexural(0);   % n = 0
assert(testHermiticity(gew.op));
gew = guide.flexural(1);   % n = 1, generates complex matrices! 
assert(testHermiticity(gew.op));

function ret = testHermiticity(op)
    h2 = ishermitian(op.L2); 
    h1 = ishermitian(op.L1, 'skew');
    h0 = ishermitian(op.L0);
    m0 = ishermitian(op.M);
    ret = h2 && h1 && h0 && m0;
end
