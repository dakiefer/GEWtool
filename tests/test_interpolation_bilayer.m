% Test whether field interpolation is consistent with the nodal values
%
% Run using: runtests(). see also: 
% https://fr.mathworks.com/help/matlab/matlab_prog/write-script-based-unit-tests.html
%
% 2024 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France

mat1 = Material('zircaloy'); % transversely isotropic
mat2 = Material('chromium'); % isotropic
h1 = 0.7e-3; h2 = 0.3e-3;   % thicknesses
N1 = 10; N2 = 8;            % number of nodes
k = 5/(h1+h2);              % compute at one wavenumber
indk = 1;                   % index of wavenumber
indw = 2;                   % we will inspect this mode
plate = Plate([mat1, mat2], [h1, h2], [N1, N2]);
gew = plate.Lamb;
dat = computeW(gew,k,indw);
datMode = extractModes(dat,indk,indw);
u = datMode.u; 
T = stress(gew, datMode); 

% % %%%% Retrieve modal field:  %%%%%%
% nodal points:
y1 = gew.geom.z{1}; y2 = gew.geom.z{2};  % nodal points of the layers
y = [y1; y2]; 

% retrieve displacements: 
u1 = squeeze(u{1}(1,1,:,:)); % displacements of layer 1
u2 = squeeze(u{2}(1,1,:,:)); % displacements of layer 2
ux = [u1(:,1); u2(:,1)];
uy = [u1(:,2); u2(:,2)];

% retrieve the y-traction is ey*T = [Tyx, Tyy]. This is continuous across the
% layer interface and is zero at the boundaries. 
T1 = squeeze(T{1}(1,1,:,2,:));
T2 = squeeze(T{2}(1,1,:,2,:));
Tyx = [T1(:,1); T2(:,1)];
Tyy = [T1(:,2); T2(:,2)];

% interpolated displacements u and stress components:
[ui, yi] = GEWinterpolate(gew, u, 1000); 
uix = ui(:,1); uiy = ui(:,2);
Ti = GEWinterpolate(gew, T, yi);
Tiyx = Ti(:,2,1); Tiyy = Ti(:,2,2);

%% test matching of interpolated and nodal values
tol = 1e-2; % tolerance for testing
ydiffs = yi - y.'; % each column is difference to nodal coordinates ''y''
[ydiffmin, indyi] = min(abs(ydiffs));
% test matching close to nodal points:
diff = ( abs(ux - uix(indyi)) )/norm(ux);
assert( all(diff < tol) )
if exist('show', 'var') && show, fprintf("ux max-difference to interpolation: %g %%\n", max(diff)*1e2); end
diff = ( abs(uy - uiy(indyi)) )/norm(uy);
assert( all(diff < tol) )
if exist('show', 'var') && show, fprintf("uy max-difference to interpolation: %g %%\n", max(diff)*1e2); end
diff = ( abs(Tyx - Tiyx(indyi)) )/norm(Tyx);
assert( all(diff < tol) )
if exist('show', 'var') && show, fprintf("Tyx max-difference to interpolation: %g %%\n", max(diff)*1e2); end
diff = ( abs(Tyy - Tiyy(indyi)) )/norm(Tyy);
assert( all(diff < tol) )
if exist('show', 'var') && show, fprintf("Tyy max-difference to interpolation: %g %%\n", max(diff)*1e2); end

%% test boundary conditions
tol = 1e-3; % tolerance for testing
tractionBottom = squeeze(T1(1,:)) / norm( T1(:) );
assert( all(abs(tractionBottom) < tol) )
tractionTop = squeeze(T2(end,:)) / norm( T2(:) );
assert( all(abs(tractionTop) < tol) )

%% test interface conditions
tol = 1e-3; % tolerance for testing
% displacement continuity:
Unorm = norm( [u1(:); u2(:)] ); 
displ1 = u1(end,:) / Unorm; 
displ2 = u2(1,:) / Unorm; 
diff = abs( displ2 - displ1 ); 
assert( all(diff < tol) )
% traction continuity:
Tnorm = norm( [T1(:); T2(:)] ); % common norm for both layers
traction1 = squeeze(T1(end,:)) / Tnorm;
traction2 = squeeze(T2(1,:)) / Tnorm;
diff = abs( traction2 - traction1 ); 
assert( all(diff < tol) )



% % %%%%%%% plot if requested %%%%%%%%
if exist('show', 'var') && show
    % All displacement components need to be continuous across the layer interfaces.
    figure(5), clf
    subplot(2,1,1); hold on, % plot ux
    plot(real(uix), yi/1e-3, 'k-');
    plot(imag(uix), yi/1e-3, 'r-');
    plot(real(ux),  y/1e-3, 'k*');
    plot(imag(ux),  y/1e-3, 'r*');
    yline(gew.geom.zItf(2)/1e-3, '-', {mat2.name, mat1.name}, 'LineWidth', 1,...
        'LabelVerticalAlignment', 'middle')
    ylim([gew.geom.zItf(1), gew.geom.zItf(end)]/1e-3)
    xlabel('modal displacement ux'), ylabel('y in mm')
    legend({'real ux', 'imag ux'}, 'Location','best')
    
    subplot(2,1,2); hold on % plot uy
    plot(real(uiy), yi/1e-3, 'k-');
    plot(imag(uiy), yi/1e-3, 'r-');
    plot(real(uy),  y/1e-3, 'k*');
    plot(imag(uy),  y/1e-3, 'r*');
    yline(gew.geom.zItf(2)/1e-3, '-', {mat2.name, mat1.name}, 'LineWidth', 1,...
        'LabelVerticalAlignment', 'middle')
    ylim([gew.geom.zItf(1), gew.geom.zItf(end)]/1e-3)
    xlabel('modal displacement uy'), ylabel('y in mm')
    legend({'real uy', 'imag uy'}, 'Location','best')
    
    % plot also the modal stresses
    % The stress components Tyx and Tyy need to be continuous across the layers. Moreover, 
    % they vanish at the boundaries. 
    figure(6), clf,
    subplot(2,1,1); hold on, % plot Tyx
    plot(real(Tiyx), yi/1e-3, 'k-');
    plot(imag(Tiyx), yi/1e-3, 'r-');
    plot(real(Tyx), y/1e-3, 'k*');
    plot(imag(Tyx), y/1e-3, 'r*');
    yline(gew.geom.zItf(2)/1e-3, '-', {mat2.name, mat1.name}, 'LineWidth', 1,...
        'LabelVerticalAlignment', 'middle')
    ylim([gew.geom.zItf(1), gew.geom.zItf(end)]/1e-3)
    xlabel('modal stress Tyx'), ylabel('y in mm')
    legend({'real Tyx', 'imag Tyx'}, 'Location','best')
    
    subplot(2,1,2); hold on, % plot Tyy
    plot(real(Tiyy), yi/1e-3, 'k-');
    plot(imag(Tiyy), yi/1e-3, 'r-');
    plot(real(Tyy), y/1e-3, 'k*');
    plot(imag(Tyy), y/1e-3, 'r*');
    yline(gew.geom.zItf(2)/1e-3, '-', {mat2.name, mat1.name}, 'LineWidth', 1,...
        'LabelVerticalAlignment', 'middle')
    ylim([gew.geom.zItf(1), gew.geom.zItf(end)]/1e-3)
    xlabel('modal stress Tyy'), ylabel('y in mm')
    legend({'real Tyy', 'imag Tyy'}, 'Location','best')
end