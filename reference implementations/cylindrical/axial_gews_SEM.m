%% Axially guided waves in a cylindrical tube 
% This script demonstrats the idea of how waves in a cylinder are computed using
% the Spectral Element Method (SEM). We represent and integrate the basis
% polynomials using chebfun, see https://www.chebfun.org 
% 
% Note: this implementation is not optimized for efficiency. The integration of
% the element matrices is much faster in the actual GEWtool implementation. 
% 
% see also:
% [1] A. Marzani, E. Viola, I. Bartoli, F. Lanza di Scalea, and P. Rizzo, “A
% semi-analytical finite element formulation for modeling stress wave
% propagation in axisymmetric damped waveguides,” Journal of Sound and
% Vibration, vol. 318, no. 3, pp. 488–505, 2008, doi: 10.1016/j.jsv.2008.04.028.
% [2] M. Zheng, C. He, Y. Lyu, and B. Wu, “Guided waves propagation in anisotropic 
% hollow cylinders by Legendre polynomial solution based on state-vector formalism,” 
% Composite Structures, vol. 207, pp. 645–657, Jan. 2019, 
% doi: 10.1016/j.compstruct.2018.09.042.
% [3] J. Zemanek Jr., “An Experimental and Theoretical Investigation of Elastic
% Wave Propagation in a Cylinder,” The Journal of the Acoustical Society of
% America, vol. 51, no. 1B, pp. 265–283, Jan. 1972, doi: 10.1121/1.1912838.
% 
% 2022-2024 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France

% specify parameters:
a = 5e-3;       % inner radius
h = 1e-3;       % thickness of the shell
n = 0;          % flexural order (circumferential wavenumber)
N = 15;         % number of collocation points
k = linspace(1e-2, 15, 70).'/h; % wavenumber-thickness (solve for frequency)
rho = 7900; lbd = 1.1538e11; mu = 7.6923e10; % material parameters

% verify that chebfun is installed:
if ~exist('chebfun','file')
    error('You need to have chebfun installed to use this example. See https://www.chebfun.org');
end
if a == 0
    error('You cannot use this code to compute solid cylinders but you can set the inner radius very small.');
end

% % derived parameters and normalization:
b = a + h;  % outer radius
II = eye(3).*shiftdim(eye(3), -2); % 4th order "unit tensor"
c = lbd*II + mu*(permute(II, [1 3 4 2]) + permute(II, [1 3 2 4])); % stiffness tensor
c0 = c(1,2,1,2); h0 = 1.7e-3; % normalization parameters
rho0 = rho; f0 = sqrt(c0/rho0)/h0; % normalization parameters
rhon = rho/rho0; cn = c/c0; hn = h/h0; % normalized
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
% polynomial finite element basis represented by chebfuns:
[yi, ~] = lobpts(N, [-1, 1]); % BUG in chebfun: does only work on dom = [-1 1]!
yi = yi/2 + 1/2;              % scale to [-1/2, 1/2]
Psi = chebfun.lagrange(yi);   % Lagrange polynomials on Gauss-Lobatto points
Psid = diff(Psi);             % differentiated polynomials
% % ALTERNATIVE: use chebyshev polynomials for FEM-expansion:
% % [yi, wi] = chebpts(2*N, [0 1]);
% % Psi = chebpoly(0:N-1, [0 1]); % Chebyshev polynomials on Chebyshev points
% % Psid = diff(Psi);             % differentiated polynomials
% compute "element matrices" (numerical integration): 
eta = chebfun('x', [0, 1]); % computational coordinate eta in [0, 1].
rnfun = a/h + eta;          % normalized radial coordinates: [r1 r2]/h
obj.PP = elemPP(Psi);
obj.PPd = elemPPd(Psi, Psid);
obj.PPr = elemPPr(Psi, rnfun);
obj.PPdr = elemPPdr(Psi, Psid, rnfun);
obj.PdPdr = elemPdPdr(Psid, rnfun);
obj.PPinvr = elemPPinvr(Psi, rnfun);

%% problem setup: (i*kh)^2*L2 + (i*kh)*L1(i*n) + L0(i*n) + w^2*M = 0
% element stiffness:
K2 = kron(cxx, obj.PPr)*hn^2;
K1 = kron(cxr, obj.PPdr)*hn + kron( cxp*(1i*n*I + A) + (1i*n*I + A)*cpx , obj.PP)*hn;
K0 = kron((1i*n*I + A)*cpr, obj.PPd) + kron((1i*n*I + A)*cpp*(1i*n*I + A), obj.PPinvr);
% element flux:
G1 = kron( crx , -obj.PPdr.' )*hn;
G0 = kron( crr , -obj.PdPdr.' )  +  kron( crp*(1i*n*I + A) , -obj.PPd.' );
% combine to polynomial of (ik):
L2 = K2; L1 = K1 + G1; L0 = K0 + G0;
M = kron(rhon*I, obj.PPr)*hn^2;

%% solve for frequency and plot:
whn = nan(length(k), size(M, 2)); tic 
khn = k*h0; % re-scale to computational units
for ii = 1:length(khn)
    wh2 = eig(khn(ii)^2*L2 - 1i*khn(ii)*L1 - L0, M); % faster if second arg is positive definite
    whn(ii,:) = sort(sqrt(wh2));
end
w = real(whn)*f0; % w in SI-Units
chron = toc; fprintf('nF: %d, nK: %d, elapsed time: %g, time per point: %g. ms\n', size(w, 1), size(w, 2), chron, chron/length(w(:))*1e3);

% % reference to test against: 
mat = MaterialIsotropic('noname', lbd, mu, rho);
guide = Cylinder(mat, [a, b], N);
gew = guide.fullyCoupled(n); 
dat = computeW(gew, k); 

% plot:
figure(1); clf; hold on;
phref = plot(dat.k, dat.w/2/pi, 'x', 'Color', 0.7*[1 1 1]); 
kk = k.*ones(size(w)); % expand to same size
phSEM = plot(kk, w/2/pi,  'k');
ylim([0, 6e3/h]);
legend([phSEM(1), phref(1)], {'SEM', 'GEWtool'}, 'Location', 'southeast')
ylabel('frequency in Hz'), xlabel('wavenumber in rad/m')


%% functions to integrate element matrices:
function me = elemPP(P) 
    % elemPP: integral of product matrix of ansatz functions ∫P*Pdy (element mass)
    me = zeros(size(P,2));
    for i = 1:size(P,2)
        for j = i:size(P,2)
            me(i,j) = sum(P(:,i)*P(:,j));
            me(j,i) = me(i,j);
        end
    end
end
function le1 = elemPPd(P, Pd) 
    % elemPPd: integral of product matrix of ∫P*P'dy (element stiffness and flux)
    le1 = zeros(size(P,2));
    for i = 1:size(P,2)
        for j = 1:size(Pd,2)
            le1(i,j) = sum(P(:,i)*Pd(:,j));
        end
    end
end
function ppr = elemPPr(P, r) 
    % elemPPr: integral ∫P*P*r dr of basis functions P
    ppr = zeros(size(P,2));
    for i = 1:size(P,2)
        for j = i:size(P,2)
            ppr(i,j) = sum(r*P(:,i)*P(:,j));
            ppr(j,i) = ppr(i,j);
        end
    end
end
function ppdr = elemPPdr(P, Pd, r) 
    % elemPPdr: integral ∫P*P'*1/r dr of basis functions P
    ppdr = zeros(size(P,2));
    for i = 1:size(P,2)
        for j = 1:size(Pd,2)
            ppdr(i,j) = sum(r*P(:,i)*Pd(:,j));
        end
    end
end
function k0 = elemPdPdr(Pd, r) 
    % elemPdPd: integral of product matrix of ∫P'*P'dy (element flux)
    k0 = zeros(size(Pd,2));
    for i = 1:size(Pd,2)
        for j = i:size(Pd,2)
            k0(i,j) = sum(r*Pd(:,i)*Pd(:,j));
            k0(j,i) = k0(i,j);
        end
    end
end
function ppr2 = elemPPinvr(P, r) 
    % elemPPr2: integral ∫P*P*1/r^2 dr of basis functions P
    ppr2 = zeros(size(P,2));
    for i = 1:size(P,2)
        for j = i:size(P,2)
            ppr2(i,j) = sum((1/r)*P(:,i)*P(:,j));
            ppr2(j,i) = ppr2(i,j);
        end
    end
end
