% compute leaky waves

mat = MaterialIsotropic('brass');   % load from database (or create your own)
mat.E = mat.E*(1 + 1e-6i);
h = 1e-3;                        % thickness in m
N = 15;                          % number of nodes (dictates accuracy)
w = 2*pi*linspace(1e-3, 3, 200)*1e6;   % wavenumbers to solve for
plate = Plate(mat, h, N);        % create waveguide description 
gewFree = plate.Lamb; tic;        % choose S+A Lamb waves (assembles matrices)
nModes = 20; 

datFree = computeK(gewFree,w,nModes);
ceFree = energyVelAxial(datFree);
PxFree = powerFluxAxial(datFree); 
pFree = poyntingVec(datFree); pFree = pFree{1};

%% incorporate fluid loading
gew = copy(gewFree);
L2 = gew.op.L2; L1 = gew.op.L1; L0 = gew.op.L0; M = gew.op.M;
rho0 = gew.np.rho0; fh0 = gew.np.fh0;
rhof = 1000/rho0; cf = 1500/fh0; % water loading

N = gew.geom.N;
dofA = 2*N + 1; % dof of degree of freedom for the fluid A
dofB = 2*N + 2; % dof of degree of freedom for the fluid B
dofuA = N+1;      % dof of uz displacement in plate that is in contact with A
dofuB = 2*N;        % dof of uz displacement in plate that is in contact with B
nDof = 2*N + 2; % size of matrices

% expand matrices: 
L2(nDof,nDof) = 0;
L1(nDof,nDof) = 0;
L0(nDof,nDof) = 0;
M(nDof,nDof) = 0;
R1 = zeros(nDof);
R2 = zeros(nDof);

% continuity of normal displacements: ibeta*A - uy = 0
L0(dofA,dofuA) = -1; 
R1(dofA,dofA)  = +1; 
L0(dofB,dofuB) = -1;
R2(dofB,dofB)  = +1;

% balance of tractions: add the boundary term [v*tA] to the FE matrices. 
% the traction induced by the fluid is tA = -w^2*rhoA*ey*A
M(dofuA,dofA) = -rhof;
M(dofuB,dofB) = +rhof;

% % assemble state space matrices
R = R2 - R1;
Z = zeros(size(R));

LL2 = [L2  ,   -R   ;
       Z   ,   L2    ];
LL1 = blkdiag(L1,L1);
LL0 = [L0  ,   Z    ; 
       R   ,   L0   ];
MM  = [M   ,   -1/cf^2*R  ; 
       Z   ,      M       ];


% gew.op.L2 = LL2; gew.op.L1 = LL1; gew.op.L0 = LL0; gew.op.M = MM;
% clear opts; opts.eigenvecs = true; opts.standardEVP = false; opts.subspace = false;
% dat = computeK(gew, w, 30, opts); toc; % solve and save 4 modes (argument optional)
dat = solveLeaky(LL2,LL1,LL0,MM,w,2*nModes,gew.np);
dat.gew = gew;
Px = powerFluxAxial(dat);
ce = energyVelAxial(dat);
p =  poyntingVec(dat); p = p{1}; 
P0 = mean(abs(Px(:)));
% p = p./abs(Px);
p = p/max(abs(p(:,:,end,1)),[],'all');
pzTop = p(:,:,end,2);
pzBottom = p(:,:,1,2);

%% plot
% indLeaky = dat.k >= -inf;
% tol = 0;
% indLeaky = pzTop >= -tol & pzBottom <= tol;
% indLeaky = abs(imag(dat.k)) > 1e-5;
indLeaky = abs(dat.k/h) < 1e2;

nModesPowerFlux = length(find(indLeaky)); 
nModesTotal = length(find(~isnan(dat.k)));
fprintf('%d modes are leaky or trapped of %d total.\n', nModesPowerFlux, nModesTotal);

figure(1); clf; hold on
plot(real(dat.k(indLeaky))/1e3, dat.w(indLeaky)/2/pi/1e6, '.','SeriesIndex',1,'DisplayName','leaky');
xlabel('wavenumber k in rad/mm'), ylabel('frequency f in MHz')
xlim([-1 1]*8)

indk = 9; indw = 42;
figure(1); plot(real(dat.k(indk,indw))/1e3, dat.w(indk,indw)/2/pi/1e6, 'rd','DisplayName','selection');
figure(33); clf; hold on; 
y = gew.geom.y{1};
% uFree = squeeze(datFree.u{1}(indk,indw,:,:)); 
uLeaky = squeeze(dat.u{1}(indk,indw,:,:)); 
subplot(2,1,1); cla; hold on;
plot(y, [real(uLeaky(:,1)), imag(uLeaky(:,1))]); 
subplot(2,1,2); cla; hold on;
plot(y, [real(uLeaky(:,2)), imag(uLeaky(:,2))]);

figure(2); clf; hold on
plot(dat.w(indLeaky)/2/pi/1e6, imag(dat.k(indLeaky))/1e3, '.','SeriesIndex',1,'DisplayName','leaky');
xlabel('frequency f in MHz'); ylabel('Im k in rad/mm')
ylim([0, 0.3])

% figure(3); clf; hold on
% plot(datFree.w/2/pi/1e6, pFree(:,:,1,2), '.','SeriesIndex',2,'DisplayName','free');
% plot(dat.w(indLeaky)/2/pi/1e6, pyTop(indLeaky), '.','SeriesIndex',1,'DisplayName','leaky');
% xlabel('frequency f in MHz'); ylabel('pyTop'); legend(legendUnq)

% figure(5); clf; hold on
% plot(datFree.w/2/pi/1e6, PxFree, '.','SeriesIndex',2,'DisplayName','free');
% plot(dat.w(indLeaky)/2/pi/1e6, Px(indLeaky), '.','SeriesIndex',1,'DisplayName','leaky');
% xlabel('frequency f in MHz'); ylabel('Px'); legend(legendUnq)
% 
% figure(6); clf; hold on
% plot(datFree.w/2/pi/1e6, ceFree, '.','SeriesIndex',2,'DisplayName','free');
% plot(dat.w(indLeaky)/2/pi/1e6, ce(indLeaky), '.','SeriesIndex',1,'DisplayName','leaky');
% xlabel('frequency f in MHz'); ylabel('ce'); legend(legendUnq)

function dat = solveLeaky(LL2,LL1,LL0,MM,w,nModes,np)
    tic
    wn = w*np.h0/np.fh0;
    kh = nan(nModes,length(wn));
    betaPsiPsi = zeros(nModes,length(wn),size(MM,1));
    for i = 1:numel(wn)
        LL0d = LL0 + wn(i)^2*MM;
        [betaPsiPsin, ikhn] = polyeig(LL0d, LL1, LL2); % compute 1i*kh: double as fast than computing kh
        [khn,ind] = sort(-1i*ikhn); % extract kh
        betaPsiPsin = betaPsiPsin(:,ind);  % sort as khn
        kh(:,i) = khn(1:nModes);
        betaPsiPsi(:,i,:) = betaPsiPsin(:,1:nModes).';
    end
    toc
    dat.k = kh/np.h0;
    dat.w = w.*ones(size(kh));
    N = (size(MM,1)/2-2)/2;
    indSecondBlock = (2*N+2+1):(2*(2*N+2));
    Psi = betaPsiPsi(:,:,indSecondBlock);
    AB = Psi(:,:,2*N+1:2*N+2);
    uInside = Psi(:,:,1:2*N);
    dat.u{1} = reshape(uInside,nModes,length(wn),N,2);
    dat.A = AB;
end