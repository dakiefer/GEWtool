% GEWintegrate: Compare integration per layer vs. integration of all layers
%
% Run using: runtests(). see also: 
% https://fr.mathworks.com/help/matlab/matlab_prog/write-script-based-unit-tests.html
%
% 2024 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France

mat1 = Material('steel');
mat2 = Material('silicon');
k = linspace(1e-2, 12, 20)/1e-3;
plate = Plate([mat1 mat2], [0.7e-3 0.3e-3], 6);
gew = plate.Lamb;
dat = computeW(gew, k, 3);
px = poyntingVecAxial(dat); % quantity to integrate
Pref = powerFluxAxial(dat); % reference

%% test without argument 
P = GEWintegrate(gew,px);
err = (P - Pref)./Pref;
assert( max(err,[],'all') <= 1e-6 )

%% test with empty argument
P = GEWintegrate(gew,px,[]);
err = (P - Pref)./Pref;
assert( max(err,[],'all') <= 1e-6 )

%% test with explicit layer IDs
P = GEWintegrate(gew,px,1:gew.geom.nLay);
err = (P - Pref)./Pref;
assert( max(err,[],'all') <= 1e-6 )

%% test each layer separately
P1 = GEWintegrate(gew,px,1);
P2 = GEWintegrate(gew,px,2);
P = P1 + P2;
err = (P - Pref)./Pref;
assert( max(err,[],'all') <= 1e-6 )
I = GEWintegrateEachLayer(gew,px);
assert(all(I{1} == P1,'all') && all(I{2} == P2,'all'));
