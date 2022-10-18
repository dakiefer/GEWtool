%% Dispersion calculation in a generally anisotropic elastic plate
% Implements the spectral element method to compute guided elastic waves.
% Depends on chebfun (https://www.chebfun.org) to represent interpolation
% functions and perform differentiation/integration.
%
% see also:
% D. A. Kiefer, "Elastodynamic quasi-guided waves for transit-time ultrasonic flow
% metering", ser. FAU Forschungen, Reihe B, Medizin, Naturwissenschaft, Technik, 
% vol. 42. Erlangen: FAU University Press, 2022, doi: 10.25593/978-3-96147-550-6
% 
% 2022 - Daniel A. Kiefer
% Institut Langevin, ESPCI Paris, France
%

h = 1e-3;   % thickness in m
N = 10;     % discretization: polynomial order of interpolants
rho = 7900; lbd = 1.1538e11; mu = 7.6923e10; % steel material

% define and normalize parameters:
II = eye(3).*shiftdim(eye(3), -2); % 4th order "unit tensor"
c = lbd*II + mu*(permute(II, [1 3 4 2]) + permute(II, [1 3 2 4])); % stiffness tensor
c0 = c(1,2,1,2); rho0 = rho; fh0 = sqrt(c0/rho0);  % normalization parameters
rhon = rho/rho0; cn = c/c0; % normalize

% relevant material matrices: 
udof = 1:2; % Lamb and/or SH
cxx = squeeze(cn(1,udof,udof,1));
cxy = squeeze(cn(1,udof,udof,2));
cyx = squeeze(cn(2,udof,udof,1));
cyy = squeeze(cn(2,udof,udof,2));
I = eye(size(cxx)); 

% % discretize: 
% % using Lagrange polynomials on GLL points:
[yi] = lobpts(N, [-1, 1]); % does only work on dom = [-1 1]!!!!
yi = yi/2; % scale to [-1/2, 1/2]
P = chebfun.lagrange(yi);
Pd = diff(P);

% % "element" matrices for one displacement component:
pp = elemPP(P);         % ∫ Pi*Pj dy
pd = elemPPd(P, Pd);    % ∫ Pi*Pj' dy
dd = elemPdPd(Pd);      % ∫ Pi'*Pj' dy
% % assemble for the displacement components:
M  = kron(rhon*I,pp);
L2 = kron(cxx, pp);
L1 = kron(cxy, pd) - kron(cyx, pd.');
L0 = kron(cyy, dd);

%% solve for frequency and plot:
kh = linspace(1e-2, 15, 300); % wavenumber*thickness 
whn = nan(size(M, 2), length(kh)); tic 
for ii = 1:length(kh)
    kh0 = kh(ii);
    [wh2] = eig(-(1i*kh0)^2*L2 - (1i*kh0)*L1 - L0, M, 'chol');  % use Cholesky: positive definite B
    whn(:,ii) = sort(sqrt(wh2));
end
fh = whn/2/pi*fh0;
chron = toc; fprintf('nF: %d, nK: %d, elapsed time: %g, time per point: %g. ms\n', size(fh, 1), size(fh, 2), chron, chron/length(fh(:))*1e3);

% % plot wavenumbers:
kkh = kh.*ones(size(fh));
figure(1), clf, hold on, plot(kkh(:)/h/1e3, real(fh(:))/h/1e6, '.');
xlim([0, 12]), ylim([0, 6]),
xlabel('k in rad/mm'), ylabel('f in MHz'),

%% solve for wavenumbers:
fh = linspace(1e-2, 6, 300).'*1e6*h; % frequency*thickness
kh = nan(length(fh), size(M, 2)*2); tic
for ii = 1:length(fh)
    whn = 2*pi*fh(ii)/fh0; % current frequency-thickness (normalized)
    [un, khi] = polyeig(L0 + whn^2*M, L1, L2); 
    kh(ii, :) = -1i*khi;
end
kh(abs(kh) >= 14.3) = nan; % only small wavenumbers are numerically accurate
chron = toc; fprintf('nF: %d, nK: %d, elapsed time: %g, time per point: %g. ms\n', size(kh, 1), size(kh, 2), chron, chron/length(kh(:))*1e3);

% % plot wave numbers:
ffh = fh.*ones(size(kh));
figure(3), clf, hold on,
scatter3(real(kh(:))/h/1e3, imag(kh(:))/h/1e3, ffh(:)/h/1e6, 8, abs(imag(kh(:))), 'filled'), 
caxis([0, 0.01]), xlim([-12, 12]), ylim([-12, 12]), view(0,0);
xlabel('real k in rad/mm'), ylabel('imag k in rad/mm'), zlabel('f in MHz')

%% Group velocity:
% % solve for frequencies as previously but save also the eigenvectors:
kh = linspace(1e-2, 15, 300); % wavenumber*thickness 
whn = nan(size(M, 2), length(kh));
u = nan([size(whn), size(M, 1)]); tic 
for ii = 1:length(kh)
    kh0 = kh(ii); 
    [ui, wh2] = eig(-(1i*kh0)^2*L2 - (1i*kh0)*L1 - L0, M, 'chol', 'vector');  % use Cholesky: positive definite B
    whn(:,ii) = sqrt(wh2);
    u(:,ii,:) = ui.';
end
fh = whn/2/pi*fh0;
chron = toc; fprintf('nF: %d, nK: %d, elapsed time: %g, time per point: %g. ms\n', size(fh, 1), size(fh, 2), chron, chron/length(fh(:))*1e3);

% % compute group velocity:
kkh = kh.*ones(size(fh)); % expand wavenumbers to matrix of same size as fh
cgn = zeros(size(kkh));
for i = 1:size(kkh,1)
    for j = 1:size(kkh,2)
        u0 = squeeze(u(i,j,:));
        cgn(i,j) = (u0'*(2*kkh(i,j).*L2 - 1i*L1)*u0) ./ (u0'*2*whn(i,j)*M*u0);
    end
end
cg = real(cgn)*fh0;

figure(15), clf, hold on
plot(fh(:)/h/1e6, cg(:)/1e3, '.'), ylim([-2, 6]), xlim([0, 6])
xlabel('f in MHz'), ylabel('cg in mm/us')

%% element matrices:
function pp = elemPP(P) 
    % elemPP - integral ∫P*Pdy of basis functions P (element mass)
    N = size(P,2);
    pp = zeros(N);
    for i = 1:N
        for j = i:N
            pp(i,j) = sum(P(:,i)*P(:,j));
            pp(j,i) = pp(i,j);
        end
    end
end

function pd = elemPPd(P, Pd) 
    % elemPPd - integral ∫P*P'dy of basis functions P (element stiffness and flux)
    N = size(P,2);
    pd = zeros(N);
    for i = 1:N
        for j = 1:N
            pd(i,j) = sum(P(:,i)*Pd(:,j));
        end
    end
end

function dd = elemPdPd(Pd)
    % elemPdPd - integral ∫P'*P'dy of basis functions P (element flux)
    N = size(Pd,2);
    dd = zeros(N);
    for i = 1:N
        for j = 1:N
            dd(i,j) = -sum(Pd(:,i)*Pd(:,j));
            dd(j,i) = dd(i,j);
        end
    end
end

