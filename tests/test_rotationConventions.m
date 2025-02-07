% Run using: runtests()
% Test the rotation conventions, i.e., extrinsic/intrinsic, passive/active.
%
% GEWtool implements intrinsic and passive rotations in rotateEuler() and
% Material.rotateEuler(). The rotation matrices are obtained by
% eulerAnglesToRotationMatrix(), where these conventions are fixed: 
% - passive/active: sign of the angle in the rotation matrix
% - intrinsic/extrinsic: order in which successive rotations are performed
% 
% Rule to convert extrinsic to intrinsic: exchange order
% extrinsic Ry(a)*Rx(a) == intrinsic Rx(a)*Ry(a)
% 
% Rule to convert passive to active: invert total rotation matrix
% passive Ry(a)*Rx(a) == Rx(-a)*Ry(-a)
% 
% see also: 
% https://fr.mathworks.com/help/matlab/matlab_prog/write-script-based-unit-tests.html
%
% 2025 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France

v = [0;0;1]; % a test vector to rotate; 
vx_PassiveIntrinsic = [0;1;0]; % 90°x in passive-intrinsic
vxy_PassiveIntrinsic = [ 0; 1; 0]; % 90°x then 90°y in passive-intrinsic
vxy_PassiveExtrinsic = [-1; 0; 0]; % 90°x then 90°y in passive-extrinsic
vxy_ActiveIntrinsic  = [ 1; 0; 0]; % 90°x then 90°y in active-intrinsic
vxy_ActiveExtrinsic  = [ 0;-1; 0]; % 90°x then 90°y in active-extrinsic

%% is passive x
% coordinate system is turned rather than the vector itself
vnew = round( rotateEuler(v, pi/2, 'x'), 8 ); 
assert(all(vnew == vx_PassiveIntrinsic))

%% passive-intrinsic xy
% You have the passive-intrinsic sequence Ry(a)*Rx(a). 
% In passive-intrinsic this is Ry(a)*Rx(a). 
vnew = round( rotateEuler(v, pi/2, 'x', pi/2, 'y'), 8);
assert(all(vnew == vxy_PassiveIntrinsic)) % second rotation does nothing: vector || y

%% passive-extrinsic xy
% You have the passive-extrinsic sequence Ry(a)*Rx(a). 
% In passive-intrinsic this is Rx(a)*Ry(a). 
vnew = round( rotateEuler(v, pi/2, 'y', pi/2, 'x'), 8);
assert(all(vnew == vxy_PassiveExtrinsic))

%% active-intrinsic xy
% You have the active-intrinsic sequence Ry(a)*Rx(a). 
% In passive-intrinsic this is Rx(-a)*Ry(-a). 
vnew = round( rotateEuler(v, -pi/2, 'y', -pi/2, 'x'), 8);
assert(all(vnew == vxy_ActiveIntrinsic))

%% active-extrinsic xy
% You have the active-extrinsic sequence Ry(a)*Rx(a). 
% In passive-intrinsic this is Ry(-a)*Rx(-a). 
vnew = round( rotateEuler(v, -pi/2, 'x', -pi/2, 'y'), 8);
assert(all(vnew == vxy_ActiveExtrinsic))

%% options "extrinsic" | "active"
% Test X-cut, propagation in Y+30°-direction of LiNbO3: 
mat0 = MaterialPiezoelectric('lithium_niobate');
% There are four ways to rotate the material to the wave propagation system.
% Assume the initial coordinate system is xyz, which becomes x'y'z' after the
% first rotation. Then, the four rotation conventions are: 
% 1. passive and intrinsic: 90° y, then 90°+30° z':
mat1 = mat0.rotateEuler(pi/2, 'y', (90+30)/180*pi, 'z'); % passive, intrinsic per default
% 2. passive and extrinsic: 90° y, then 90°+30° x:
mat2opt = mat0.rotateEuler(pi/2, 'y', (90+30)/180*pi, 'x', "extrinsic"); % passive default
mat2 = mat0.rotateEuler((90+30)/180*pi, 'x', pi/2, 'y'); % passive, extrinsic converted to passive intrinsic
% 3. active and intrinsic: -90° y, then -90°-30° x': 
mat3opt = mat0.rotateEuler(-pi/2, 'y', -(90+30)/180*pi, 'x', "active"); % intrinsic default
mat3 = mat0.rotateEuler((90+30)/180*pi, 'x', pi/2, 'y'); % active, intrinsic converted to passive intrinsic
% 4. active and extrinsic: -90° y, then -90°-30° z: 
mat4opt = mat0.rotateEuler(-pi/2, 'y', -(90+30)/180*pi, 'z', "active", "extrinsic");
mat4 = mat0.rotateEuler(pi/2, 'y', (90+30)/180*pi, 'z'); % active, extrinsic converted to passive intrinsic

assert(mat2opt == mat1 && mat2 == mat1);
assert(mat3opt == mat1 && mat3 == mat1); 
assert(mat4opt == mat1 && mat4 == mat1); 