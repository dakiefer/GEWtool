% % EDAT:
mat = Material.Material('brass'); 
h = 1e-3; % plate thickness in mm
N = 40; % number of collocation points
plate = Dispersion.Plate(mat, h, N);

% compute wavenumbers and wave displacement structures:
freq = linspace(0.01, 4, 150).'*1e6;
Nmodes = 60;
data = plate.solveForFreq(freq, Nmodes);
ff = freq.*ones(size(data.k)); % for plotting
[ySCM, wSCM] = chebpts(N, [-0.5 0.5]*h);
data.u = normalizeL2(data.u, wSCM);  % normalize to ∫ conj(u).u dy = 1

% % new implementation:
clear plate
plate = Plate(mat, h, N);
gew = plate.Lamb;
% k = linspace(1e-3, 15, 150)/h;
dat = computeK(gew, 2*pi*freq);
[yNew, wNew] = Layer.nodes(N); wNew = wNew*h; yNew = (yNew-0.5)*h;
dat.u{1} = normalizeL2(dat.u{1}, wNew); % % normalize to ∫ conj(u).u dy = 1

%% displacements:
uOld = abs(data.u(:,:,1,1));
uNew = abs(dat.u{1}(:,:,1,1));
clf, hold on, title('displacement ux at boundary')
plot(ff(:), abs(uOld(:)), 'x')
plot(dat.w(:)/2/pi, abs(uNew(:)), '.')

%% particle velocity:
vOld = abs(data.v(:,:,1,1));
vNew = velocity(dat);
vNew = abs(vNew{1});
vNew = vNew(:,:,1,1);
figure, hold on, title('velocity vx at boundary')
plot(ff(:), real(vOld(:)), 'x')
plot(dat.w(:)/2/pi, real(vNew(:)), '.')

%% strain: 
Sold = abs(data.S(:,:,1,1,1));
Snew = strain(gew, dat);
Snew = abs(Snew{1});
Snew = Snew(:,:,1,1,1);
figure, hold on, title('strain Sxx at boundary')
plot(ff(:), real(Sold(:)), 'x')
plot(dat.w(:)/2/pi, Snew(:), '.')
% ylim([0, 5000])

%% stress: 
TOld = abs(data.T(:,:,1,1,1));
Tnew = stress(gew, dat);
Tnew = abs(Tnew{1});
Tnew = Tnew(:,:,1,1,1);
figure, hold on, title('stress Txx at boundary')
plot(ff(:), real(TOld(:)), 'x')
plot(dat.w(:)/2/pi, real(Tnew(:)), '.'), %ylim([0, 3e14])

%% power flux at interface 
pOld = data.poyntingVectors;
pOld = pOld(:,:,1,1);
pNew = poyntingVec(gew, dat);
pNew = pNew{1}(:,:,1,1);
figure, hold on, title('power flux density px at boundary')
plot(ff(:), real(pOld(:)), 'x')
plot(dat.w(:)/2/pi, real(pNew(:)), '.')

%% kinetik energy 
EkOld = data.Ek;
EkNew = energyKinetic(gew, dat);
figure, hold on, title('total kinetik energy')
plot(ff(:), real(EkOld(:)), 'x')
plot(dat.w(:)/2/pi, real(EkNew(:)), '.')

%% elastic energy:
EsOld = data.Es;
EsNew = energyElastic(gew, dat);
figure, hold on, title('total elastic energy')
plot(ff(:), abs(EsOld(:)), 'x');
plot(dat.w(:)/2/pi, abs(EsNew(:)), '.');

%% power flux total
Pold = data.Px;
Pnew = powerFlux(gew, dat);
figure, hold on, title('total power flux P')
plot(ff(:), real(Pold(:)), 'x')
plot(dat.w(:)/2/pi, real(Pnew(:)), '.')

%% energy vel
ceOld = data.ce;
ceNew = energyVel(gew, dat);
figure, hold on, title('energy velocity ce')
plot(ff(:), ceOld(:), 'x')
plot(dat.w(:)/2/pi, ceNew(:), '.')
