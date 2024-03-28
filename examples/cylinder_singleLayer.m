%% Compute axially guided waves in a tube
% Different circumferential orders n of the waves are compared.
% 
% 2022-2024 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France

% % specify parameters:
mat = MaterialIsotropic('steel');       % load from material/database
r = [10e-3, 11e-3]; h = r(end) - r(1);  % radii and thickness 
N = 8;                                  % number of nodes
k = linspace(0, 4, 200)/h; % wavenumbers to solve for
cyl = Cylinder(mat, r, N); % create waveguide description 

%% compute and plot
n = 0:1:2; % circumferential order to solve for
figure, hold on; ax = gca; xlim([0, k(end)/1e3]); clear ph;
xlabel('wavenumber k in rad/mm'), ylabel('frequency f in MHz')
cc = winter(length(n));
for i = 1:length(n)
    gew = cyl.fullyCoupled(n(i)); % assembles matrices for the specified waves
    dat = computeW(gew, k); 
    phs = plot(dat.k/1e3, dat.w/2/pi/1e6, 'Color', cc(i,:)); ph(i) = phs(1);
    ylim([0, 2]); drawnow;
end
legend(ph, strcat('n = ', num2str(n.')), 'Location', 'southeast')
