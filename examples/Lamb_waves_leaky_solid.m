% Compute leaky waves in a brass plate coupled on one side to teflon. 

mat.lambda = 81.3120e9; mat.mu = 40.6560e9; mat.rho = 8.4e3; % brass parameters (does not correspond to GEWtool parameters)
mat = MaterialIsotropic(mat);
h = 1e-3;                        % thickness in m
N = 18;                          % number of nodes (dictates accuracy)
w = 2*pi*linspace(1e-3, 3, 150)*1e6;   % wavenumbers to solve for
plate = Plate(mat, h, N);        % create waveguide description 
gewFree = plate.Lamb; tic;        % choose S+A Lamb waves (assembles matrices)
nModes = 25; 

datFree = computeK(gewFree,w,nModes);
ceFree = energyVelAxial(datFree);
pFree = poyntingVec(datFree); pFree = pFree{1};

%% incorporate fluid loading
gew = copy(gewFree);
rho0 = gew.np.rho0; fh0 = gew.np.fh0; c0 = gew.np.c0;
matA.lbd =  2.6785e9/c0; matA.mu = 0.6655e9/c0; matA.rho = 2.2e3/rho0; % teflon
II = eye(3).*shiftdim(eye(3), -2); % 4th order "unit tensor"
matA.c = matA.lbd*II + matA.mu*(permute(II, [1 3 4 2]) + permute(II, [1 3 2 4])); % stiffness tensor
matA.ct = sqrt(matA.mu/matA.rho); matA.cl = sqrt((matA.lbd + 2*matA.mu)/matA.rho); % wave speeds

[L2,L1,L0,M,Rkb,Rke,Rl,Rt] = couplingSolid(gew, matA);
Z = zeros(size(L0));
al = 1/(matA.cl^2); at = 1/(matA.ct^2);
M = M - al*Rl - at*Rt; 

% % assemble state space matrices
LL3 = [ Z,   -Rkb,  -Rke,  Z  ;
        Z,    Z,     Z,   -Rke; 
        Z,    Z,     Z,   -Rkb; 
        Z,    Z,     Z,    Z ];
LL2 = blkdiag(L2, L2, L2, L2);
LL1 = [ L1,   Z,     Z,    Z  ; 
        Rkb,  L1,    Z,    Z  ; 
        Rke   Z,     L1,   Z  ; 
        Z,    Rke,   Rkb,  L1];
LL0 = blkdiag(L0, L0, L0, L0);
MM  = blkdiag(M, M, M, M);
MM1 = [ Z,    -al*Rkb, -at*Rke,   Z   ;
        Z,     Z,      Z,      -at*Rke;
        Z,     Z,      Z,      -al*Rkb; 
        Z,     Z,      Z,         Z  ]; 


dat = solveLeaky(LL3,LL2,LL1,LL0,MM,MM1,w,4*nModes,gew.np);
Px = powerFluxAxial(dat);
ce = energyVelAxial(dat);
p =  poyntingVec(dat); p = p{1}; 
p = p/max(abs(p(:,:,1,2)),[],'all');
% pyTop = p(:,:,end,2);
pyBottom = p(:,:,1,2);

%% validate
% addpath('~/Projekte/Radiation/leaky wave calc')
% folder = fileparts(which('getMatricesFEMsolid2'));
% fileName = fullfile(folder, 'results','oldResultsSBFEM','plate_BrassTeflon','compare_cp_fewer_points.fig');
% open(fileName); hold on;
% % indRef = abs(imag(kRef)) < 0.5*abs(real(kRef)); % & real(betah) > 1e-4 & real(etah) > 1e-4;
% % plot(ffRef(indRef), 2*pi*ffRef(indRef)./real(kRef(indRef)), 'ob', 'MarkerSize',6,'DisplayName','multipar');
% ind = imag(dat.k*h) > 1e-3 & imag(dat.k*h) < 0.5*abs(real(dat.k*h)); % & real(betah) < -1e-4 & real(etah) < -1e-4; % is empty!
% plot(dat.w(ind)/2/pi/1e6, dat.w(ind)./real(dat.k(ind))/1e3, 'xg', 'MarkerSize',4,'LineWidth',1,'DisplayName','state');
% xlim([0, w(end)/2/pi/1e6])
% legend; drawnow;
% % % plot attenuation
% fileName = fullfile(folder, 'results','oldResultsSBFEM','plate_BrassTeflon','compare_att_fewer_points.fig');
% open(fileName);
% hold all
% % att=imag(kRef)/h*20/log(10)*1000;
% % plot(ffRef(indRef), att(indRef), 'ob', 'MarkerSize',6,'DisplayName','MultiPar');
% att=imag(dat.k)*20/log(10);
% plot(dat.w(ind)/2/pi/1e6, att(ind), 'xg', 'MarkerSize',4,'LineWidth',1,'DisplayName','state');
% xlim([0, w(end)/2/pi/1e6])
% drawnow;

%% plot
% indLeaky = dat.k >= -inf;
tol = 1e-6;
indLeaky = pyBottom <= tol;
% indLeaky = abs(imag(dat.k)) < 1e-1;

nModesPowerFlux = length(find(indLeaky)); 
nModesTotal = length(find(~isnan(dat.k)));
fprintf('%d modes are leaky or trapped of %d total.\n', nModesPowerFlux, nModesTotal);

figure(1); clf; hold on
plot(real(dat.k)/1e3, dat.w/2/pi/1e6, '.','SeriesIndex',2,'DisplayName','all');
plot(real(dat.k(indLeaky))/1e3, dat.w(indLeaky)/2/pi/1e6, '.','SeriesIndex',1,'DisplayName','leaky');
xlabel('wavenumber k in rad/mm'), ylabel('frequency f in MHz')
xlim([-1 1]*8)

figure(76); clf; hold on
plot3(real(dat.k)/1e3, imag(dat.k)/1e3, dat.w/2/pi/1e6, '.','SeriesIndex',2,'DisplayName','all');
plot3(real(dat.k(indLeaky))/1e3, imag(dat.k(indLeaky))/1e3, dat.w(indLeaky)/2/pi/1e6, '.','SeriesIndex',1,'DisplayName','leaky');
xlabel('Re k'), ylabel('Im k'), zlabel('f'); legend(legendUnq)
xlim([-1 1]*8), ylim([-1 1]*1.5)

indk = 9; indw = 42;
figure(1); plot(real(dat.k(indk,indw))/1e3, dat.w(indk,indw)/2/pi/1e6, 'rd','DisplayName','selection');
figure(33); clf; hold on; 
y = gew.geom.y{1};
uLeaky = squeeze(dat.u{1}(indk,indw,:,:)); 
subplot(2,1,1); cla; hold on;
plot(y, [real(uLeaky(:,1)), imag(uLeaky(:,1))]); 
subplot(2,1,2); cla; hold on;
plot(y, [real(uLeaky(:,2)), imag(uLeaky(:,2))]);

figure(2); clf; hold on
plot(dat.w(indLeaky)/2/pi/1e6, imag(dat.k(indLeaky))/1e3, '.','SeriesIndex',1,'DisplayName','leaky');
xlabel('frequency f in MHz'); ylabel('Im k in rad/mm')
ylim([0, 0.3])

figure(3); clf; hold on
plot(dat.w/2/pi/1e6, pyBottom, '.','SeriesIndex',2,'DisplayName','all');
plot(dat.w(indLeaky)/2/pi/1e6, pyBottom(indLeaky), '.','SeriesIndex',1,'DisplayName','leaky');
xlabel('frequency f in MHz'); ylabel('pyBottom'); legend(legendUnq)

figure(5); clf; hold on
plot(dat.w(indLeaky)/2/pi/1e6, Px(indLeaky), '.','SeriesIndex',1,'DisplayName','leaky');
xlabel('frequency f in MHz'); ylabel('Px'); legend(legendUnq)

figure(6); clf; hold on
plot(datFree.w/2/pi/1e6, ceFree, '.','SeriesIndex',2,'DisplayName','free');
plot(dat.w(indLeaky)/2/pi/1e6, ce(indLeaky), '.','SeriesIndex',1,'DisplayName','leaky');
xlabel('frequency f in MHz'); ylabel('ce'); legend(legendUnq)



function dat = solveLeaky(LL3,LL2,LL1,LL0,MM,MM1,w,nModes,np)
    tic
    NPsi = size(LL0,1)/4;
    forthBlock = 3*NPsi+1:4*NPsi;
    wn = w*np.h0/np.fh0;
    indModes = 8+(1:nModes); % the first couple of eigenvalues are zero!
    kh = nan(nModes,length(wn));
    Psi = zeros(nModes,length(wn),NPsi);
    for i = 1:numel(wn)
        % kappal = wn(i)/(matA.cl); % only needed to compute transverse wavenumbers
        % kappat = wn(i)/(matA.ct); % only needed to compute transverse wavenumbers
        LL0d = LL0 + wn(i)^2*MM;
        LL1d = LL1 + wn(i)^2*MM1;
        [eVeci, ikhn] = polyeig(LL0d, LL1d, LL2, LL3); % compute 1i*kh: double as fast than computing kh
        [khn,ind] = sort(-1i*ikhn); % extract kh
        eVeci = eVeci(:,ind);  % sort as khn
        kh(:,i) = khn(indModes);
        Psi(:,i,:) = eVeci(forthBlock,indModes).';
    end
    toc
    dat.k = kh/np.h0;
    dat.w = w.*ones(size(kh));
    N = (size(MM,1)/4-2)/2;
    AB = Psi(:,:,2*N+1:2*N+2);
    uInside = Psi(:,:,1:2*N);
    dat.u{1} = reshape(uInside,nModes,length(wn),N,2);
    dat.A = AB;
end



function [L2,L1,L0,M,R1,R2,R3,R4] = couplingSolid(gew, matUnb)
    % Code adapted from Hauke's code getMatricesFEM() lines 812-987
    
    L2 = gew.op.L2; L1 = gew.op.L1; L0 = gew.op.L0; M = gew.op.M;
    udof = 1:2;
    N = gew.geom.N;
    dofA = 2*N + (1:2); % degree of freedom for longitudinal and transverse bulk wave in exterior
    dofuA = [1, N+1]; % dofs of displacements in plate that are in contact with A
    nDof = 2*N + 2; % new size of matrices
    
    % expand matrices: 
    L2(nDof,nDof) = 0;
    L1(nDof,nDof) = 0;
    L0(nDof,nDof) = 0;
    M(nDof,nDof) = 0;
    
    % initialized coupling matrices
    % coupling is of the form
    % k^2*R0 + k*alpha*r1 + k*beta*r2 + kl^2*r3 + ks^2*r4
    R0 = zeros(nDof);
    R1 = zeros(nDof);
    R2 = zeros(nDof);
    R3 = zeros(nDof);
    R4 = zeros(nDof);
    
    % get material properties
    cyx = squeeze(matUnb.c(2,udof,udof,1));
    cyy = squeeze(matUnb.c(2,udof,udof,2));

    % helper matrices
    D1 = diag([1 0]);
    D2 = diag([0 1]);
    A0 = [1 0; 0 -1];
    A1 = [0 0; 1 0];
    A2 = [0 1; 0 0];
    
    r01 = -(cyx*A0 - cyy*(A1*D1+A2*D2));
    r11 = -(cyx*A1 + cyy*A0*D1);
    r21 = -(cyx*A2 + cyy*A0*D2);
    r31 = -cyy*A1*D1;
    r41 = -cyy*A2*D2;
    
    % contributions to unbounded equation, bottom
    signA = -1;
    L1(dofA,dofuA) = -eye(2);
    R0(dofA,dofA) =  A0;
    R1(dofA,dofA) = -A1;
    R2(dofA,dofA) = -A2;
    R0(dofuA,dofA) = -signA*r01;
    R1(dofuA,dofA) =  signA*r11;
    R2(dofuA,dofA) =  signA*r21;
    R3(dofuA,dofA) =  signA*r31;
    R4(dofuA,dofA) =  signA*r41;
    
    % change sign to adjust to Daniel's notation, including a factor i^2
    R1 = -R1;
    R2 = -R2;
    R3 = -R3;
    R4 = -R4;
    
    L2 = L2 + R0; % coupling to k^2 term
end
