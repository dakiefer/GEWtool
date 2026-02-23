% Test that the energy velocity is correctly computed in piezoelectric waveguides. 
% 
% This compares the axial energy velocity cex to the axial group velocity cgx. A
% comparison to dw/dk is also done. 
%
% Run using: runtests(). see also: 
% https://fr.mathworks.com/help/matlab/matlab_prog/write-script-based-unit-tests.html
%
% 2026 - Daniel A. Kiefer, Institut Langevin, CNRS, ESPCI Paris, France

h = 1e-3;
N = 10;
k = linspace(1e-2, 15, 100)/h;
tol = 1e-10;

%% AlN - Lamb: energy velocity == group velocity
mat = MaterialPiezoelectric('aluminum_nitrate');
% matYX = mat.rotateEuler(-pi/2,'x'); % passive intrinsic rotation
guide = Plate(mat, h, N);
gew = guide.Lamb; 
gew = gew.fixGdof(gew.geom.gdofBC{1}(3,:)); % fix potential at top and bottom -> closed boundary conditions! 
dat = computeW(gew, k, 10);

cgx = groupVelAxial(dat);
cex = energyVelAxial(dat);

dev = abs(cgx - cex)/mat.cl;
if exist('show', 'var') && show
    maxDev = max(dev,[],'all')
    figure(3); clf; hold on; title('AlN - Lamb')
    plot(dat.w(:), cgx(:), 'rx', 'DisplayName','cg');
    plot(dat.w(:), cex(:), 'k.', 'DisplayName','ce');
    legend('Location','southeast');
end
assert( max(dev,[],'all') <= tol )

%% AlN - fullyCoupled: energy velocity == group velocity
mat = MaterialPiezoelectric('aluminum_nitrate');
mat = mat.rotateEuler(-pi/7,'x'); % some strange orientation
guide = Plate(mat, h, N);
gew = guide.fullyCoupled; 
gew = gew.fixGdof(gew.geom.gdofBC{1}(4,:)); % fix potential at top and bottom -> closed boundary conditions! 
dat = computeW(gew, k, 10);

cgx = groupVelAxial(dat);
cex = energyVelAxial(dat);

dev = abs(cgx - cex)/mat.cl;
if exist('show', 'var') && show
    maxDev = max(dev,[],'all')
    figure(3); clf; hold on; title('AlN - Lamb')
    plot(dat.w(:), cgx(:), 'rx', 'DisplayName','cg');
    plot(dat.w(:), cex(:), 'k.', 'DisplayName','ce');
    legend('Location','southeast');
end
assert( max(dev,[],'all') <= tol )

%% AlN - fullyCoupled: dwdk == group velocity
mat = MaterialPiezoelectric('lithium_niobate');
mat = mat.rotateEuler(-pi/5,'x'); % some strange orientation (S-A are coupled)
mat = mat.rotateEuler(-pi/7,'z'); % some strange orientation (Lamb-sh are coupled)
guide = Plate(mat, h, N);
gew = guide.fullyCoupled; 
gew = gew.fixGdof(gew.geom.gdofBC{1}(4,:)); % fix potential at top and bottom -> closed boundary conditions! 
k = linspace(1e-2, 10, 1000)/h;
dat = computeW(gew, k, 7);

cgx = groupVelAxial(dat);
dwdk = diff(dat.w)./diff(dat.k); dwdk(end+1,:) = nan;

dev = abs(cgx - dwdk)/mat.cl;
if exist('show', 'var') && show
    maxDev = max(dev,[],'all')
    figure(3); clf; hold on; title('LNB - fullyCoupled')
    plot(dat.w(:), cgx(:), 'rx', 'DisplayName','cg');
    plot(dat.w(:), dwdk(:), 'k.', 'DisplayName','dwdk');
    legend('Location','southeast');
end
assert( max(dev,[],'all') <= 2e-2 ) % less strick checking: dwdk is not accurate

%% LNB - fullyCoupled: energy velocity == group velocity
mat = MaterialPiezoelectric('lithium_niobate');
mat = mat.rotateEuler(-pi/7,'x'); % some strange orientation
guide = Plate(mat, h, N);
gew = guide.fullyCoupled; 
gew = gew.fixGdof(gew.geom.gdofBC{1}(4,:)); % fix potential at top and bottom -> closed boundary conditions! 
dat = computeW(gew, k, 10);

cgx = groupVelAxial(dat);
cex = energyVelAxial(dat);

dev = abs(cgx - cex)/mat.cl;
if exist('show', 'var') && show
    maxDev = max(dev,[],'all')
    figure(3); clf; hold on; title('fully coupled')
    plot(dat.w(:), cgx(:), 'rx', 'DisplayName','cg');
    plot(dat.w(:), cex(:), 'k.', 'DisplayName','ce');
    legend('Location','southeast');
end
assert( max(dev,[],'all') <= tol )
