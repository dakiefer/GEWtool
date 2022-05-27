%% EDAT:
mat = Material.Material('brass'); 
h = 1e-3; % plate thickness in mm
N = 25; % number of collocation points
plate = Dispersion.Plate(mat, h, N);

% compute wavenumbers and wave displacement structures:
freq = linspace(0.01, 4, 150).'*1e6;
Nmodes = 60;
data = plate.solveForFreq(freq, Nmodes);
ff = freq.*ones(size(data.k)); % for plotting

%% new implementation:
clear plate
plate = Plate(mat, h, N);
wguide = plate.Lamb;
% k = linspace(1e-3, 15, 150)/h;
dat = computeK(wguide, 2*pi*freq);

%% displacements:
uOld = abs(data.u(:,:,1,1));
uNew = abs(dat.u{1}(:,:,1,1));
figure, hold on, title('displacement ux at boundary')
plot(ff(:), real(uOld(:)), 'x')
plot(dat.w(:)/2/pi, real(uNew(:)), '.')

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
Snew = strain(wguide, dat);
Snew = abs(Snew{1});
Snew = Snew(:,:,1,1,1);
figure, hold on, title('strain Sxx at boundary')
plot(ff(:), real(Sold(:)), 'x')
plot(dat.w(:)/2/pi, Snew(:), '.')
ylim([0, 5000])

%% stress: 
TOld = abs(data.T(:,:,1,1,1));
Tnew = stress(wguide, dat);
Tnew = abs(Tnew{1});
Tnew = Tnew(:,:,1,1,1);
figure, hold on, title('stress Txx at boundary')
plot(ff(:), real(TOld(:)), 'x')
plot(dat.w(:)/2/pi, real(Tnew(:)), '.'), ylim([0, 3e14])

%% kinetik energy 
EkOld = data.Ek;
EkNew = energyKinetic(wguide, dat);
figure, hold on, title('total kinetik energy')
plot(ff(:), real(EkOld(:)), 'x')
plot(dat.w(:)/2/pi, real(EkNew(:)), '.')

%% elastic energy:
EsOld = data.Es;
EsNew = energyElastic(wguide, dat);
figure, hold on, title('total elastic energy')
plot(ff(:), abs(EsOld(:)), 'x');
plot(dat.w(:)/2/pi, abs(EsNew(:)), '.');

%% power flux at interface 
pOld = data.poyntingVectors;
pOld = pOld(:,:,1,1);
pNew = poyntingVec(wguide, dat);
pNew = pNew{1}(:,:,1,1);
figure, hold on, title('power flux density px at boundary')
plot(ff(:), real(pOld(:)), 'x')
plot(dat.w(:)/2/pi, real(pNew(:)), '.')

%% power flux total
Pold = data.Px;
Pnew = powerFlux(wguide, dat);
figure, hold on, title('total power flux P')
plot(ff(:), real(Pold(:)), 'x')
plot(dat.w(:)/2/pi, real(Pnew(:)), '.')

%% energy vel
ceOld = data.ce;
ceNew = energyVel(wguide, dat);
figure, hold on, title('energy velocity ce')
plot(ff(:), ceOld(:), 'x')
plot(dat.w(:)/2/pi, ceNew(:), '.')
