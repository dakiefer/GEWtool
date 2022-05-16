%% Dispersion calculation in a generally anisotropic elastic plate
% 
% 2022 - Daniel A. Kiefer
% Institut Langevin, Paris, France
% 

h = 1e-3;   % thickness in m
N = 19;     % discretization: polynomial order of interpolants
rho = 7900; lbd = 1.1538e11; mu = 7.6923e10; % steel material

% define and normalize parameters:
II = eye(3).*shiftdim(eye(3), -2); % 4th order "unit tensor"
c = lbd*II + mu*(permute(II, [1 3 4 2]) + permute(II, [1 3 2 4])); % stiffness tensor
c0 = c(1,2,1,2); h0 = h; % normalization parameters
rho0 = rho; fh0 = sqrt(c0/rho0); % normalization parameters
rhon = rho/rho0; cn = c/c0;

% relevant material matrices: 
udof = 1:3; % Lamb and/or SH
cxx = squeeze(cn(1,udof,udof,1));
cxy = squeeze(cn(1,udof,udof,2));
cyx = squeeze(cn(2,udof,udof,1));
cyy = squeeze(cn(2,udof,udof,2));
I = eye(size(cxx)); 

%% discretize: 
tic 
dom = [0 1];
[y, w] = chebpts(N); % lobpts, legpts
% D1 = diffmat(N,1,[0 1],'chebkind2'); D2 = diffmat(N,2,[0 1],'chebkind2');
% mesh = Geometry(dom, N, length(udof));
% mesh.y{1} = y;
Psi = chebpoly(0:N, dom);
Psid = diff(Psi);
me = elemM(Psi);
k2 = me;
k1 = elemK1(Psi, Psid);
g1 = -k1.';
g0 = elemG0(Psid);
M  = kron(rhon*I,me);
K2 = kron(cxx, k2);
K1 = kron(cxy, k1);
G1 = kron(cyx, g1);
G0 = kron(cyy, g0);

L2 = K2; L1 = K1 + G1; L0 = G0;
chron = toc; fprintf('assembling time: %g\n', chron);

%% solve for frequency and plot:
kh = linspace(1e-2, 15, 300); % wavenumber*thickness 
whn = nan(size(M, 2), length(kh)); tic 
for ii = 1:length(kh)
    kh0 = kh(ii);
    [wh2] = eig((1i*kh0)^2*L2 + (1i*kh0)*L1 + L0, -M, 'chol'); 
    whn(:,ii) = sqrt(wh2);
end
fh = whn/2/pi*fh0; fh(fh == 0) = nan;
chron = toc; fprintf('nF: %d, nK: %d, elapsed time: %g, time per point: %g. ms\n', size(fh, 1), size(fh, 2), chron, chron/length(fh(:))*1e3);

% plot wavenumbers:
kkh = kh.*ones(size(fh));
figure(1), hold on, plot(kkh(:)/h/1e3, real(fh(:))/h/1e6, '.');
xlim([0, 12]), ylim([0, 6]),
xlabel('k in rad/mm'), ylabel('f in MHz'),

ind = find(imag(fh)~=0);
figure(1), hold on, plot(kkh(ind)/h/1e3, real(fh(ind))/h/1e6, 'x');
drawnow

%% solve for wavenumbers:
fh = linspace(1e-2, 6, 300).'*1e6*h; % frequency*thickness
kh = nan(length(fh), size(M, 2)*2); tic
for ii = 1:length(fh)
    whn = 2*pi*fh(ii)/fh0; % current frequency-thickness (normalized)
    [un, khi] = polyeig(L0 + whn^2*M, L1, L2); 
    kh(ii, :) = -1i*khi;
end
chron = toc; fprintf('nF: %d, nK: %d, elapsed time: %g, time per point: %g. ms\n', size(kh, 1), size(kh, 2), chron, chron/length(kh(:))*1e3);

% % plot wave numbers:
ffh = fh.*ones(size(kh));
figure(2), scatter(real(kh(:))/h/1e3, ffh(:)/h/1e6, 8, abs(imag(kh(:))), 'filled'), 
caxis([0, 0.12]), xlim([0, 12]), ylim([0, fh(end)/h/1e6])
xlabel('k in rad/mm'), ylabel('f in MHz')


%% element matrices:
function me = elemM(Psi) 
    me = zeros(length(Psi));
    for i = 1:length(Psi)
        for j = i:length(Psi)
            me(i,j) = sum(Psi(:,i)*Psi(:,j));
            me(j,i) = me(i,j);
        end
    end
end

function le1 = elemK1(Psi, Psid) 
    le1 = zeros(size(Psi,2));
    for i = 1:size(Psi,2)
        for j = 1:size(Psi,2)
            le1(i,j) = sum(Psi(:,i)*Psid(:,j));
        end
    end
end

function g0 = elemG0(Psid) 
    g0 = zeros(size(Psid,2));
    for i = 1:size(Psid,2)
        for j = i:size(Psid,2)
            g0(i,j) = -sum(Psid(:,i)*Psid(:,j));
            g0(j,i) = g0(i,j);
        end
    end
end

