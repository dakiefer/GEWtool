%% Axially guided waves in a cylindrical tube 
% The spectral collocation method (SCM) is implemented to solve the boundary-value
% problem describing waves guided along the axis of a cylindrical tube.
%
% See also 
% [1] A. T. I. Adamou and R. V. Craster, “Spectral methods for modelling guided
% waves in elastic media,” The Journal of the Acoustical Society of America,
% vol. 116, no. 3, pp. 1524–1535, Sep. 2004, doi: 10.1121/1.1777871.
% [2] D. A. Kiefer, Elastodynamic quasi-guided waves for transit-time ultrasonic
% flow metering. Erlangen: FAU University Press, 2022. 
% doi: 10.25593/978-3-96147-550-6.
% [3] J. Zemanek Jr., “An Experimental and Theoretical Investigation of Elastic
% Wave Propagation in a Cylinder,” The Journal of the Acoustical Society of
% America, vol. 51, no. 1B, pp. 265–283, Jan. 1972, doi: 10.1121/1.1912838.
% 
% 2022-2024 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France
%             Pierre Chantelot, Institut Langevin, ESPCI Paris, France

% specify parameters:
a = 5e-3;       % inner radius
h = 1e-3;       % thickness of the shell
n = 0;          % flexural order (circumferential wavenumber)
N = 15;         % number of collocation points
k = linspace(1e-2, 15, 70).'/h; % wavenumber-thickness (solve for frequency)
rho = 7900; lbd = 1.1538e11; mu = 7.6923e10; % material parameters

% verify that chebfun is installed:
if ~exist('chebdif','file')
    error('You need to have DMSUITE installed to use this example. See https://fr.mathworks.com/matlabcentral/fileexchange/29-dmsuite');
end
if a == 0
    error('You cannot use this code to compute solid cylinders but you can set the inner radius very small.');
end

% % derived parameters and normalization:
b = a + h; % outer radius
II = eye(3).*shiftdim(eye(3), -2); % 4th order "unit tensor"
c = lbd*II + mu*(permute(II, [1 3 4 2]) + permute(II, [1 3 2 4])); % stiffness tensor
c0 = c(1,2,1,2); rho0 = rho; fh0 = sqrt(c0/rho0); % normalization parameters
rhon = rho/rho0; cn = c/c0; % normalized
% relevant material matrices: 
cxx = squeeze(cn(1,:,:,1));
crr = squeeze(cn(2,:,:,2));
cpp = squeeze(cn(3,:,:,3));
crp = squeeze(cn(2,:,:,3));
cpr = squeeze(cn(3,:,:,2));
cxp = squeeze(cn(1,:,:,3));
cpx = squeeze(cn(3,:,:,1));
crx = squeeze(cn(2,:,:,1));
cxr = squeeze(cn(1,:,:,2));
A = [0, 0, 0; 0, 0, -1; 0, 1, 0]; % differetiation in curvilinear coordinate system
I = eye(size(A));

%% discretization 
[y_dash, Dr_dash] = chebdif(N, 2);  % y_dash in [-1, 1]
r = (h*y_dash + b + a)/2;           % radial coordinate
rn = r/h;                           % uses thickness h as normalization thickness
Dr1 =  2*Dr_dash(:,:,1);            % differentiation on unit domain
Dr2 = (2)^2*Dr_dash(:,:,2);         % 2nd order differentiation on unit domain
Id = eye(size(Dr1));                % identity matrix for discretization
rn1inv = diag(1./rn);               % necessary for BC

%% problem setup: (i*kh)^2*L2 + (i*kh)*L1(i*n) + L0(i*n) + w^2*M = 0
L2 = kron(cxx, Id); 
L1 = kron(cxr + crx, Dr1) + kron(crx + cxp*(1i*n*I + A) + (1i*n*I + A)*cpx, diag(1./rn));
L0 =   kron( crr, Dr2 ) ...
     + kron( crr + crp*(1i*n*I + A) + (1i*n*I + A)*cpr, diag(1./rn)*Dr1 ) ...
     + kron( (1i*n*I + A)*cpp*(1i*n*I + A), diag(1./rn.^2) ); 
M = kron(rhon*I, Id);

% incorporate BCs:
B1 = kron(crx,  Id([1, N], :)); % BC going into L1
B0 = kron(crr, Dr1([1, N], :)) + kron(crp*(1i*n*I + A), rn1inv([1, N], :)); % into L0
dofBC = [1, N, N+1, 2*N, 2*N+1, 3*N]; % equations to replace by BCs
L2(dofBC, :) = 0; L1(dofBC, :) = B1; L0(dofBC, :) = B0; M(dofBC, :) = 0;

%% solve for frequency:
whn = nan(length(k),size(M, 2)); tic,
khn = k*h; % re-scale to computational units (use h as unit distance)
for ii = 1:length(khn)
    [wh2] = eig(khn(ii)^2*L2 - 1i*khn(ii)*L1 - L0, M, 'vector');
    whn(ii,:) = sort(sqrt(wh2));
end
w = real(whn)*fh0/h; % convert back to SI units
w(w == 0) = nan;  % remove zero eigenvalues due to singular matrix M
chron = toc; fprintf('nF: %d, nK: %d, elapsed time: %g, time per point: %g. ms\n', size(w, 2), size(w, 1), chron, chron/length(w(:))*1e3);

% % reference to test against: 
mat = MaterialIsotropic('noname', lbd, mu, rho);
guide = Cylinder(mat, [a, b], N);
gew = guide.fullyCoupled(n); 
dat = computeW(gew,k); 

% plot:
figure(1); clf; hold on;
phref = plot(dat.k, dat.w/2/pi, 'x', 'Color', 0.7*[1 1 1]); 
kk = k.*ones(size(w)); % expand to same size
phSCM = plot(kk, w/2/pi,  'k');
ylim([0, 6e3/h]);
legend([phSCM(1), phref(1)], {'SCM', 'GEWtool'}, 'Location', 'southeast')
ylabel('frequency in Hz'), xlabel('wavenumber in rad/m')