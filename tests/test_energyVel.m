% Test that the energy velocity is correctly computed by energyVelAxial(). 
% 
% This compares the axial energy velocity cex to the axial group velocity cgx.
% Both should be identical for nondissipative waveguides.
%
% Run using: runtests(). see also: 
% https://fr.mathworks.com/help/matlab/matlab_prog/write-script-based-unit-tests.html
%
% 2024 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France

h = 1e-3;
N = 10;
k = linspace(1e-2, 15, 100)/h;

%% test symmetric Lamb waves
mat = Material('steel');
guide = Plate(mat, h, N);
gew = guide.LambS;
dat = computeW(gew, k, 6); 
cgx = real(groupVelAxial(gew,dat));
cex = real(energyVelAxial(gew,dat));

devCgCe = abs(cgx - cex)/mat.cl;
assert( max(devCgCe,[],'all') <= 1e-4 )

if exist('show', 'var') && show
    maxDev = max(devCgCe,[],'all')
    figure(1); clf; hold on; title('symmetric Lamb waves')
    plot(dat.w(:), cgx(:), 'rx', 'DisplayName','cg');
    plot(dat.w(:), cex(:), 'k.', 'DisplayName','ce');
    legend('Location','southeast');
end

%% test antisymmetric Lamb waves
mat = Material('steel');
guide = Plate(mat, h, N);
gew = guide.LambA;
dat = computeW(gew, k, 6); 
cgx = real(groupVelAxial(gew,dat));
cex = real(energyVelAxial(gew,dat));

devCgCe = abs(cgx - cex)/mat.cl;
assert( max(devCgCe,[],'all') <= 1e-4 )

if exist('show', 'var') && show
    maxDev = max(devCgCe,[],'all')
    figure(2); clf; hold on; title('antisymmetric Lamb waves')
    plot(dat.w(:), cgx(:), 'rx', 'DisplayName','cg');
    plot(dat.w(:), cex(:), 'k.', 'DisplayName','ce');
    legend('Location','southeast');
end

%% test Lamb waves
mat = Material('steel');
guide = Plate(mat, h, N);
gew = guide.Lamb;
dat = computeW(gew, k, 6); 
cgx = real(groupVelAxial(gew,dat));
cex = real(energyVelAxial(gew,dat));

% S and A modes might cross. Wrong sorting leads to large errors in dwdk. Test
% only the remaining points:
devCgCe = abs(cgx - cex)/mat.cl;
assert( max(devCgCe,[],'all') <= 1e-4 )

if exist('show', 'var') && show
    maxDev = max(devCgCe,[],'all')
    figure(3); clf; hold on; title('Lamb waves, S and A coupled')
    plot(dat.w(:), cgx(:), 'rx', 'DisplayName','cg');
    plot(dat.w(:), cex(:), 'k.', 'DisplayName','ce');
    legend('Location','southeast');
end

%% test anisotropic
mat = Material('triclinic'); mat = mat.rotateEuler(0, pi/7, 0);
guide = Plate(mat, h, N);
gew = guide.fullyCoupled;
dat = computeW(gew, k, 6); 
cgx = real(groupVelAxial(gew,dat));
cex = real(energyVelAxial(gew,dat));

devCgCe = abs(cgx - cex)/mat.cl;
assert( max(devCgCe,[],'all') <= 1e-5 )

if exist('show', 'var') && show
    maxDev = max(devCgCe,[],'all')
    figure(4); clf; hold on; title('triclinic plate')
    plot(dat.w(:), cgx(:), 'rx', 'DisplayName','cg');
    plot(dat.w(:), cex(:), 'k.', 'DisplayName','ce');
    legend('Location','southeast');
end

%% test bilayer
mat = Material('steel'); 
mat2 = Material('aluminum');
guide = Plate([mat mat2], h/2, round(N/2));
gew = guide.Lamb;
dat = computeW(gew, k, 6); 
cgx = real(groupVelAxial(gew,dat));
cex = real(energyVelAxial(gew,dat));

devCgCe = abs(cgx - cex)/mat.cl;
assert( max(devCgCe,[],'all') <= 1e-4 )

if exist('show', 'var') && show
    maxDev = max(devCgCe,[],'all')
    figure(5); clf; hold on; title('Lamb waves - bilayered plate')
    plot(dat.w(:), cgx(:), 'rx', 'DisplayName','cg');
    plot(dat.w(:), cex(:), 'k.', 'DisplayName','ce');
    legend('Location','southeast');
end

%%%%%%%% NOTE %%%%%%%%%
% % Cylinder is not supported for now % %
%%%%%%%%%%%%%%%%%%%%%%%

%% test cylinder longitudinal 
mat = Material('steel');
guide = Cylinder(mat, [h, 2*h]-0.8*h/2, N); % small inner radius -> curvature is important
gew = guide.longitudinal;
dat = computeW(gew, k, 6); 
cgx = real(groupVelAxial(gew,dat));
cex = real(energyVelAxial(gew,dat));

devCgCe = abs(cgx - cex)/mat.cl;
assert( max(devCgCe,[],'all') <= 1e-4 )

% plot if the variable "show" has been set to true
if exist('show', 'var') && show
    maxDev = max(devCgCe,[],'all')
    figure(6); clf; hold on; title('Cylinder isotropic longitudinal')
    plot(dat.w(:), cgx(:), 'rx', 'DisplayName','cg');
    plot(dat.w(:), cex(:), 'k.', 'DisplayName','ce');
    legend('Location','southeast');
end

%% test cylinder n = 0
mat = Material('triclinic'); mat = mat.rotateEuler(0, pi/7, 0);
guide = Cylinder(mat, [h, 2*h]-h/2, N); % small inner radius -> curvature is important
gew = guide.fullyCoupled(0);
dat = computeW(gew, k, 6); 
cgx = real(groupVelAxial(gew,dat));
cex = real(energyVelAxial(gew,dat));

devCgCe = abs(cgx - cex)/mat.cl;
assert( max(devCgCe,[],'all') <= 1e-4 )

% plot if the variable "show" has been set to true
if exist('show', 'var') && show
    maxDev = max(devCgCe,[],'all')
    figure(7); clf; hold on; title('Cylinder triclinic n = 0')
    plot(dat.w(:), cgx(:), 'rx', 'DisplayName','cg');
    plot(dat.w(:), cex(:), 'k.', 'DisplayName','ce');
    legend('Location','southeast');
end

%% test cylinder n = 1
mat = Material('triclinic'); mat = mat.rotateEuler(0, pi/7, 0);
guide = Cylinder(mat, [h, 2*h]-h/2, N); % small inner radius -> curvature is important
gew = guide.fullyCoupled(1);
dat = computeW(gew, k, 6); 
cgx = real(groupVelAxial(gew,dat));
cex = real(energyVelAxial(gew,dat));

devCgCe = abs(cgx - cex)/mat.cl;
assert( max(devCgCe,[],'all') <= 1e-4 )

% plot if the variable "show" has been set to true
if exist('show', 'var') && show
    maxDev = max(devCgCe,[],'all')
    figure(8); clf; hold on; title('Cylinder triclinic n = 1')
    plot(dat.w(:), cgx(:), 'rx', 'DisplayName','cg');
    plot(dat.w(:), cex(:), 'k.', 'DisplayName','ce');
    legend('Location','southeast');
end
