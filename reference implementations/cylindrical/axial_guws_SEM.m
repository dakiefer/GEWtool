%% Dispersion calculation in a generally anisotropic elastic plate
% 
% 2022 - Daniel A. Kiefer
% Institut Langevin, Paris, France
% 

clear all
% reference data:
load('../../tests/data/disperse.mat')
figure(1), clf, hold on
plot(disperse.k/1e-3, disperse.f*1e6, 'g.')
xlim([0, 0.5]), ylim([0, 200])
xlabel('wavenumber k in rad/m'), ylabel('frequency f in Hz')
ax = gca; ax.ColorOrderIndex = 1; % reset color order 

% parameters:
r = [20 21];
h = r(end)-r(1); 
rho=7932; mu=rho*3260^2; lbd=rho*5960^2-2*mu;	
mat = MaterialIsotropic('noname', lbd, mu, rho);
N = 10;
kh = linspace(1e-3, 0.5, 100); % wavenumber-thickness (solve for frequency)
n = 0;

% normalized parameters
c = mat.c; rho = mat.rho;
c0 = c(1,2,1,2); h0 = h; % normalization parameters
rho0 = rho; fh0 = sqrt(c0/rho0); % normalization parameters
rhon = rho/rho0; cn = c/c0;
rn = r/h0;
hn = h/h0;
an = r(1)/h0;
detJ = @(eta) hn*(hn*eta + an);

% create matrices:
% [yi, ~] = lobpts(N, [-1, 1]); % does only work on dom = [-1 1]!!!!
% yi = yi/2 + 1/2; % scale to [-1/2, 1/2]
% Psi = chebfun.lagrange(yi);

% [yi, wi] = chebpts(2*N, [0 1]);
Psi = chebpoly(0:N-1, [0 1]); 

Psid = diff(Psi);
eta = chebfun('x', [0, 1]);
rnfun = h0*eta + an;
detJfun = chebfun(detJ, [0, 1]);
obj.PP = elemPP(Psi, detJfun);
obj.PPd = elemPPd(Psi, Psid, detJfun);
obj.PdPd = elemPdPd(Psid, detJfun);
obj.PPr = elemPPr(Psi, rnfun, detJfun);
obj.PPr2 = elemPPr2(Psi, rnfun, detJfun);
obj.PPdr = elemPPdr(Psi, Psid, rnfun, detJfun);

% relevant material matrices: 
udof = 1:3;
cxx = squeeze(cn(1,udof,udof,1));
cpp = squeeze(cn(3,udof,udof,3));
crp = squeeze(cn(2,udof,udof,3));
cpr = squeeze(cn(3,udof,udof,2));
cxp = squeeze(cn(1,udof,udof,3));
cpx = squeeze(cn(3,udof,udof,1));
cxr = squeeze(cn(1,udof,udof,2));
% differetiation in curvilinear coordinate system:
A = [0, 0, 0; 0, 0, -1; 0, 1, 0]; A = squeeze(A(udof, udof));
B = [0, 0, 0; 0, 1,  0; 0, 0, 1]; B = squeeze(B(udof, udof));
I = eye(size(A));

% element stiffness:
% K2 = kron(cxx, obj.PP);
% K1 = kron(cxr, obj.PPd) + kron( (cxp + cpx)*(1i*n*I + A) , obj.PPr);
% k0PPdr = cpp + cpr*(1i*n*I + A); % first term in K0
% k0PPr2 = -cpp*B - (crp + cpr)*A + (1i*n)*(cpp*2*A - (crp + cpr)*I) + (1i*n)^2*cpp; % second term in K0
% K0 = kron(k0PPdr, obj.PPdr) + kron(k0PPr2, obj.PPr2);
K2 = kron(cxx, obj.PP);
K1 = kron(cxr, obj.PPd) + kron( (cxp + cpx)*A , obj.PPr);
k0PPdr = cpp + cpr*A; % first term in K0
k0PPr2 = -cpp*B - (cpr + crp)*A; % second term in K0
K0 = kron(k0PPdr, obj.PPdr) + kron(k0PPr2, obj.PPr2);
% element flux:
crx = squeeze(cn(2,udof,udof,1));
crr = squeeze(cn(2,udof,udof,2));
crp = squeeze(cn(2,udof,udof,3));
% normalized element flux:
% G1 = kron( crx , -obj.PPd.' );
% G0 = kron( crr , -obj.PdPd.' )  +  kron( crp*(1i*n*I + A) , -obj.PPdr.' );
G1 = kron( crx , -obj.PPd.' );
G0 = kron( crr , -obj.PdPd.' )  +  kron( crp*A , -obj.PPdr.' );
% combine to polynomial of (ik):
L2 = K2; L1 = K1 + G1; L0 = K0 + G0;
M = kron(rhon*I, obj.PP);


%% solve for frequency and plot:
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
figure(1), hold on, plot(kkh(:)/h, real(fh(:))/h, '.');



function me = elemPP(P, detJ) 
    % elemPP: integral of product matrix of ansatz functions ∫P*Pdy (element mass)
    me = zeros(size(P,2));
    for i = 1:size(P,2)
        for j = i:size(P,2)
            me(i,j) = sum(P(:,i)*P(:,j)*detJ);
            me(j,i) = me(i,j);
        end
    end
end
function le1 = elemPPd(P, Pd, detJ) 
    % elemPPd: integral of product matrix of ∫P*P'dy (element stiffness and flux)
    le1 = zeros(size(P,2));
    for i = 1:size(P,2)
        for j = 1:size(Pd,2)
            le1(i,j) = sum(P(:,i)*Pd(:,j)*detJ);
        end
    end
end
function k0 = elemPdPd(Pd, detJ) 
    % elemPdPd: integral of product matrix of ∫P'*P'dy (element flux)
    k0 = zeros(size(Pd,2));
    for i = 1:size(Pd,2)
        for j = i:size(Pd,2)
            k0(i,j) = sum(Pd(:,i)*Pd(:,j)*detJ);
            k0(j,i) = k0(i,j);
        end
    end
end
function ppr = elemPPr(P, r, detJ) 
    % elemPPr: integral ∫P*P*1/r dr of basis functions P
    ppr = zeros(size(P,2));
    for i = 1:size(P,2)
        for j = i:size(P,2)
            ppr(i,j) = sum((1/r)*P(:,i)*P(:,j)*detJ);
            ppr(j,i) = ppr(i,j);
        end
    end
end
function ppr2 = elemPPr2(P, r, detJ) 
    % elemPPr2: integral ∫P*P*1/r^2 dr of basis functions P
    ppr2 = zeros(size(P,2));
    for i = 1:size(P,2)
        for j = i:size(P,2)
            ppr2(i,j) = sum((1/r^2)*P(:,i)*P(:,j)*detJ);
            ppr2(j,i) = ppr2(i,j);
        end
    end
end
function ppdr = elemPPdr(P, Pd, r, detJ) 
    % elemPPdr: integral ∫P*P'*1/r dr of basis functions P
    ppdr = zeros(size(P,2));
    for i = 1:size(P,2)
        for j = 1:size(Pd,2)
            ppdr(i,j) = sum((1/r)*P(:,i)*Pd(:,j)*detJ);
        end
    end
end

