%% Dispersion calculation in a generally anisotropic elastic plate
% Implements the spectral element method to solve for guided waves.
% Depends on chebfun (https://www.chebfun.org) to represent interpolation
% functions and perform differentiation/integration.
%
% 2022 - Daniel A. Kiefer
% Institut Langevin, Paris, France
% 

h = 1e-3;   % thickness in m
N = 10;     % discretization: polynomial order of interpolants
rho = 7900; lbd = 1.1538e11; mu = 7.6923e10; % steel material

% define and normalize parameters:
II = eye(3).*shiftdim(eye(3), -2); % 4th order "unit tensor"
c = lbd*II + mu*(permute(II, [1 3 4 2]) + permute(II, [1 3 2 4])); % stiffness tensor
c0 = c(1,2,1,2); h0 = h; % normalization parameters
rho0 = rho; fh0 = sqrt(c0/rho0); % normalization parameters
rhon = rho/rho0; cn = c/c0;

% relevant material matrices: 
udof = 1:2; % Lamb and/or SH
cxx = squeeze(cn(1,udof,udof,1));
cxy = squeeze(cn(1,udof,udof,2));
cyx = squeeze(cn(2,udof,udof,1));
cyy = squeeze(cn(2,udof,udof,2));
I = eye(size(cxx)); 

%% discretize: 
% % using Lagrange polynomials on GLL points:
[yi] = lobpts(N, [-1, 1]); % does only work on dom = [-1 1]!!!!
yi = yi/2; % scale to [-1/2, 1/2]
Psi = chebfun.lagrange(yi);
Psid = diff(Psi);

% % "element" matrices for one displacement component:
me = elemM(Psi);
k2 = me;
k1 = elemK1(Psi, Psid);
g1 = -k1.';
g0 = -elemG0(Psid);
% % assemble for the displacement components:
M  = kron(rhon*I,me);
K2 = kron(cxx, k2);
K1 = kron(cxy, k1);
G1 = kron(cyx, g1);
G0 = kron(cyy, g0);

L2 = K2; L1 = K1 + G1; L0 = G0;

%% solve for frequency and plot:
kh = linspace(1e-2, 15, 200); % wavenumber*thickness 
whn = nan(size(M, 2), length(kh)); tic 
for ii = 1:length(kh)
    kh0 = kh(ii);
    [wh2] = eig((1i*kh0)^2*L2 + (1i*kh0)*L1 + L0, -M); 
    whn(:,ii) = sqrt(wh2);
end
fh = whn/2/pi*fh0;
chron = toc; fprintf('nF: %d, nK: %d, elapsed time: %g, time per point: %g. ms\n', size(fh, 1), size(fh, 2), chron, chron/length(fh(:))*1e3);

% % plot wavenumbers:
kkh = kh.*ones(size(fh));
figure(1), hold on, plot(kkh(:)/h/1e3, real(fh(:))/h/1e6, '.');
xlim([0, 12]), ylim([0, 6]),
xlabel('k in rad/mm'), ylabel('f in MHz'),

%% solve for wavenumbers:
fh = linspace(1e-2, 6, 200).'*1e6*h; % frequency*thickness
kh = nan(length(fh), size(M, 2)*2); tic
for ii = 1:length(fh)
    whn = 2*pi*fh(ii)/fh0; % current frequency-thickness (normalized)
    [un, khi] = polyeig(L0 + whn^2*M, L1, L2); 
    kh(ii, :) = -1i*khi;
end
kh(abs(kh) >= 12) = nan; % only small wavenumbers are numerically accurate
chron = toc; fprintf('nF: %d, nK: %d, elapsed time: %g, time per point: %g. ms\n', size(kh, 1), size(kh, 2), chron, chron/length(kh(:))*1e3);

% % plot wave numbers:
ffh = fh.*ones(size(kh));
figure(2), hold on, scatter(real(kh(:))/h/1e3, ffh(:)/h/1e6, 8, abs(imag(kh(:))), 'filled'), 
caxis([0, 0.12]), xlim([0, 12]), ylim([0, fh(end)/h/1e6])
xlabel('k in rad/mm'), ylabel('f in MHz')
% plot3(real(kh(:))/h/1e3, imag(kh(:))/h/1e3, ffh(:)/h/1e6, 'k.'), xlim([-12, 12]), ylim([-12, 12])
% xlabel('real k in rad/mm'), ylabel('imag k in rad/mm'), zlabel('f in MHz')


%% element matrices:
function me = elemM(P) 
    N = size(P,2);
    me = zeros(N);
    for i = 1:N
        for j = i:N
            me(i,j) = sum(P(:,i)*P(:,j));
            me(j,i) = me(i,j);
        end
    end
end

function k1 = elemK1(P, Pd) 
    N = size(P,2);
    k1 = zeros(N);
    for i = 1:N
        for j = 1:N
            k1(i,j) = sum(P(:,i)*Pd(:,j));
        end
    end
end

function g0 = elemG0(Pd) 
    N = size(Pd,2);
    g0 = zeros(N);
    for i = 1:N
        for j = 1:N
            g0(i,j) = sum(Pd(:,i)*Pd(:,j));
            g0(j,i) = g0(i,j);
        end
    end
end

