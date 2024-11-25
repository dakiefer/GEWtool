% Test that the group velocity is correctly computed by groupVelAxial()
% This compares the group velocity cg computed based on the eigenvectors to the 
% finite difference approximation ∆w/∆k.
%
% Run using: runtests(). see also: 
% https://fr.mathworks.com/help/matlab/matlab_prog/write-script-based-unit-tests.html
%
% 2024 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France

h = 1e-3; 
N = 15; 
k = linspace(1e-2, 15, 1000)/h;
nModes = 6; 

%% test symmetric Lamb waves
mat = Material('steel');
guide = Plate(mat, h, N);
gew = guide.LambS;
dat = computeW(gew, k, nModes); 
cg = real(groupVelAxial(gew,dat));
dwdk = diff(dat.w,1,1)./diff(dat.k,1,1); dwdk(end+1,:) = nan; % dwdk = circshift(dwdk,1,1);

devCgDwdk = abs(cg - dwdk);
devCgDwdk = sort(devCgDwdk(:),'descend');
devCgDwdk = devCgDwdk(200:end); % remove nan entries and the biggest mismatch (large deviation due to approximation in dwdk)

% plot if the variable "show" has been set to true
if exist('show', 'var') && show
    maxDev = max(devCgDwdk,[],'all')
    figure(1); clf; hold on; title('symmetric Lamb waves')
    plot(dat.w(:), dwdk(:), 'rx', 'DisplayName','$\partial \omega / \partial k$');
    plot(dat.w(:), cg(:), 'k.', 'DisplayName','cg');
    legend('Location','southeast');
end

assert( max(devCgDwdk,[],'all') <= 50 );

%% test antisymmetric Lamb waves
mat = Material('steel');
guide = Plate(mat, h, N);
gew = guide.LambA;
dat = computeW(gew, k, nModes); 
cg = real(groupVelAxial(gew,dat));
dwdk = diff(dat.w,1,1)./diff(dat.k,1,1); dwdk(end+1,:) = nan; % dwdk = circshift(dwdk,1,1);

devCgDwdk = abs(cg - dwdk);
devCgDwdk = sort(devCgDwdk(:),'descend');
devCgDwdk = devCgDwdk(200:end); % remove nan entries and the biggest mismatch (large deviation due to approximation in dwdk)

% plot if the variable "show" has been set to true
if exist('show', 'var') && show
    maxDev = max(devCgDwdk,[],'all')
    figure(2); clf; hold on; title('antisymmetric Lamb waves')
    plot(dat.w(:), dwdk(:), 'rx', 'DisplayName','$\partial \omega / \partial k$');
    plot(dat.w(:), cg(:), 'k.', 'DisplayName','cg');
    legend('Location','southeast');
end

assert( max(devCgDwdk,[],'all') <= 50 );

%% test Lamb waves
mat = Material('steel');
guide = Plate(mat, h, N);
gew = guide.Lamb;
dat = computeW(gew, k, nModes); 
cg = real(groupVelAxial(gew,dat));
dwdk = diff(dat.w,1,1)./diff(dat.k,1,1); dwdk(end+1,:) = nan; % dwdk = circshift(dwdk,1,1);

% S and A modes might cross. Wrong sorting leads to large errors in dwdk. Test
% only the remaining points:
devCgDwdk = abs(cg - dwdk);
devCgDwdk = sort(devCgDwdk(:),'descend');
devCgDwdk = devCgDwdk(200:end); % remove nan entries and the biggest mismatch (large deviation due to approximation in dwdk)

% plot if the variable "show" has been set to true
if exist('show', 'var') && show
    maxDev = max(devCgDwdk,[],'all')
    figure(3); clf; hold on; title('Lamb waves, S and A coupled')
    plot(dat.w(:), dwdk(:), 'rx', 'DisplayName','$\partial \omega / \partial k$');
    plot(dat.w(:), cg(:), 'k.', 'DisplayName','cg');
    legend('Location','southeast');
end

assert( max(devCgDwdk,[],'all') <= 50 );

%% test anisotropic
mat = Material('triclinic'); mat = mat.rotateEuler(0, pi/7, 0);
guide = Plate(mat, h, N);
gew = guide.fullyCoupled;
dat = computeW(gew, k, nModes); 
cg = real(groupVelAxial(gew,dat));
dwdk = diff(dat.w,1,1)./diff(dat.k,1,1); dwdk(end+1,:) = nan; % dwdk = circshift(dwdk,1,1);

devCgDwdk = abs(cg - dwdk);
devCgDwdk = sort(devCgDwdk(:),'descend');
devCgDwdk = devCgDwdk(200:end); % remove nan entries and the biggest mismatch (large deviation due to approximation in dwdk)

% plot if the variable "show" has been set to true
if exist('show', 'var') && show
    maxDev = max(devCgDwdk,[],'all')
    figure(4); clf; hold on; title('triclinic plate')
    plot(dat.w(:), dwdk(:), 'rx', 'DisplayName','$\partial \omega / \partial k$');
    plot(dat.w(:), cg(:), 'k.', 'DisplayName','cg');
    legend('Location','southeast');
end

assert( max(devCgDwdk,[],'all') <= 50 );

%% test bilayer
mat = Material('steel'); 
mat2 = Material('aluminum');
guide = Plate([mat mat2], h/2, round(N/2));
gew = guide.Lamb;
dat = computeW(gew, k, nModes); 
cg = real(groupVelAxial(gew,dat));
dwdk = diff(dat.w,1,1)./diff(dat.k,1,1); dwdk(end+1,:) = nan; % dwdk = circshift(dwdk,1,1);

devCgDwdk = abs(cg - dwdk);
devCgDwdk = sort(devCgDwdk(:),'descend');
devCgDwdk = devCgDwdk(200:end); % remove nan entries and the biggest mismatch (large deviation due to approximation in dwdk)

% plot if the variable "show" has been set to true
if exist('show', 'var') && show
    maxDev = max(devCgDwdk,[],'all')
    figure(5); clf; hold on; title('Lamb waves - bilayered plate')
    plot(dat.w(:), dwdk(:), 'rx', 'DisplayName','$\partial \omega / \partial k$');
    plot(dat.w(:), cg(:), 'k.', 'DisplayName','cg');
    legend('Location','southeast');
end

assert( max(devCgDwdk,[],'all') <= 50 );

%% test cylinder longitudinal 
mat = Material('steel');
guide = Cylinder(mat, [h, 2*h]-h/2, N); % small inner radius -> curvature is important
gew = guide.longitudinal;
dat = computeW(gew, k, nModes); 
cg = real(groupVelAxial(gew,dat));
dwdk = diff(dat.w,1,1)./diff(dat.k,1,1); dwdk(end+1,:) = nan; % dwdk = circshift(dwdk,1,1);

devCgDwdk = abs(cg - dwdk);
devCgDwdk = sort(devCgDwdk(:),'descend');
devCgDwdk = devCgDwdk(200:end); % remove nan entries and the biggest mismatch (large deviation due to approximation in dwdk)

% plot if the variable "show" has been set to true
if exist('show', 'var') && show
    maxDev = max(devCgDwdk,[],'all')
    figure(6); clf; hold on; title('Cylinder isotropic longitudinal')
    plot(dat.w(:), dwdk(:), 'rx', 'DisplayName','$\partial \omega / \partial k$');
    plot(dat.w(:), cg(:), 'k.', 'DisplayName','cg');
    legend('Location','southeast');
end

assert( max(devCgDwdk,[],'all') <= 50 );

%% test cylinder n = 0
mat = Material('triclinic'); mat = mat.rotateEuler(0, pi/7, 0);
guide = Cylinder(mat, [h, 2*h]-h/2, N); % small inner radius -> curvature is important
gew = guide.fullyCoupled(0);
dat = computeW(gew, k, nModes); 
cg = real(groupVelAxial(gew,dat));
dwdk = diff(dat.w,1,1)./diff(dat.k,1,1); dwdk(end+1,:) = nan; % dwdk = circshift(dwdk,1,1);

devCgDwdk = abs(cg - dwdk);
devCgDwdk = sort(devCgDwdk(:),'descend');
devCgDwdk = devCgDwdk(200:end); % remove nan entries and the biggest mismatch (large deviation due to approximation in dwdk)

% plot if the variable "show" has been set to true
if exist('show', 'var') && show
    maxDev = max(devCgDwdk,[],'all')
    figure(7); clf; hold on; title('Cylinder triclinic n = 0')
    plot(dat.w(:), dwdk(:), 'rx', 'DisplayName','$\partial \omega / \partial k$');
    plot(dat.w(:), cg(:), 'k.', 'DisplayName','cg');
    legend('Location','southeast');
end

assert( max(devCgDwdk,[],'all') <= 50 );

%% test cylinder N = 1
mat = Material('triclinic'); mat = mat.rotateEuler(0, pi/7, 0);
guide = Cylinder(mat, [h, 2*h]-h/2, N); % small inner radius -> curvature is important
gew = guide.fullyCoupled(1);
dat = computeW(gew, k, nModes); 
cg = real(groupVelAxial(gew,dat));
dwdk = diff(dat.w,1,1)./diff(dat.k,1,1); dwdk(end+1,:) = nan; % dwdk = circshift(dwdk,1,1);

devCgDwdk = abs(cg - dwdk);
devCgDwdk = sort(devCgDwdk(:),'descend');
devCgDwdk = devCgDwdk(200:end); % remove nan entries and the biggest mismatch (large deviation due to approximation in dwdk)

% plot if the variable "show" has been set to true
if exist('show', 'var') && show
    maxDev = max(devCgDwdk,[],'all')
    figure(8); clf; hold on; title('Cylinder triclinic n = 1')
    plot(dat.w(:), dwdk(:), 'rx', 'DisplayName','$\partial \omega / \partial k$');
    plot(dat.w(:), cg(:), 'k.', 'DisplayName','cg');
    legend('Location','southeast');
end

assert( max(devCgDwdk,[],'all') <= 50 );

%% test circumferential Lamb
mat = Material('steel');
guide = CylinderCircumferential(mat, [h, 2*h]-h/2, N); % small inner radius -> curvature is important
gew = guide.Lamb;
dat = computeW(gew, k, nModes); 
cg = real(groupVelAxial(gew,dat));
dwdk = diff(dat.w,1,1)./diff(dat.k,1,1); dwdk(end+1,:) = nan; % dwdk = circshift(dwdk,1,1);

devCgDwdk = abs(cg - dwdk);
devCgDwdk = sort(devCgDwdk(:),'descend');
devCgDwdk = devCgDwdk(200:end); % remove nan entries and the biggest mismatch (large deviation due to approximation in dwdk)

% plot if the variable "show" has been set to true
if exist('show', 'var') && show
    maxDev = max(devCgDwdk,[],'all')
    figure(8); clf; hold on; title('Cylinder triclinic n = 1')
    plot(dat.w(:), dwdk(:), 'rx', 'DisplayName','$\partial \omega / \partial k$');
    plot(dat.w(:), cg(:), 'k.', 'DisplayName','cg');
    legend('Location','southeast');
end

assert( max(devCgDwdk,[],'all') <= 50 );

%% test circumferential anisotropic
mat = Material('triclinic'); mat = mat.rotateEuler(0, pi/7, 0);
guide = CylinderCircumferential(mat, [h, 2*h]-h/2, N); % small inner radius -> curvature is important
gew = guide.fullyCoupled;
dat = computeW(gew, k, nModes); 
cg = real(groupVelAxial(gew,dat));
dwdk = diff(dat.w,1,1)./diff(dat.k,1,1); dwdk(end+1,:) = nan; % dwdk = circshift(dwdk,1,1);

devCgDwdk = abs(cg - dwdk);
devCgDwdk = sort(devCgDwdk(:),'descend');
devCgDwdk = devCgDwdk(200:end); % remove nan entries and the biggest mismatch (large deviation due to approximation in dwdk)

% plot if the variable "show" has been set to true
if exist('show', 'var') && show
    maxDev = max(devCgDwdk,[],'all')
    figure(8); clf; hold on; title('Cylinder triclinic n = 1')
    plot(dat.w(:), dwdk(:), 'rx', 'DisplayName','$\partial \omega / \partial k$');
    plot(dat.w(:), cg(:), 'k.', 'DisplayName','cg');
    legend('Location','southeast');
end

assert( max(devCgDwdk,[],'all') <= 50 );