% % load pre-computed fields:
load('data/Lamb_ref_fields.mat')

% % current implementation:
N = 35;
plate = Plate(mat, h, N); % mat and h from Lamb_ref_fields.mat
gew = plate.Lamb;
dat = computeK(gew, 2*pi*freq, Nmodes); % freq and Nmodes from Lamb_ref_fields.mat
[y, wInt] = Layer.nodes(N); wInt = wInt*h; y = (y-0.5)*h;
dat.u{1} = normalizeL2(dat.u{1}, wInt); % % normalize to ∫ conj(u).u dy = 1
relTol = 1e-4;
absTol = 1e-3;

%% displacements:
u = abs(dat.u{1}(:,:,1,1));
err = (abs(u(:) - uRef(:)))./uRef(:);
assert( all(err < relTol) )

figure, hold on, title('displacement ux at boundary')
plot(datRef.w(:)/2/pi, abs(uRef(:)), 'x')
plot(dat.w(:)/2/pi, abs(u(:)), '.')

%% particle velocity:
v= velocity(dat);
v = abs(v{1});
v = v(:,:,1,1);
err = (abs(v(:) - vRef(:)))./vRef(:);
assert( all(err < relTol) )

figure, hold on, title('velocity vx at boundary')
plot(datRef.w(:)/2/pi, real(vRef(:)), 'x')
plot(dat.w(:)/2/pi, real(v(:)), '.')

%% strain: 
S = strain(gew, dat);
S = abs(S{1});
S = S(:,:,1,1,1);
err = (abs(S(:) - SRef(:)))./SRef(:);
assert( all(err < relTol) )

figure, hold on, title('strain Sxx at boundary')
plot(datRef.w(:)/2/pi, real(SRef(:)), 'x')
plot(dat.w(:)/2/pi, S(:), '.')
% ylim([0, 5000])

%% stress: 
T = stress(gew, dat);
T = abs(T{1});
T = T(:,:,1,1,1);
err = (abs(T(:) - TRef(:)))./TRef(:);
assert( all(err < relTol) )

figure, hold on, title('stress Txx at boundary')
plot(datRef.w(:)/2/pi, real(TRef(:)), 'x')
plot(dat.w(:)/2/pi, real(T(:)), '.'), %ylim([0, 3e14])

%% power flux at interface 
p = poyntingVec(gew, dat);
p = p{1}(:,:,1,1);
err = (abs(p(:) - pRef(:))); % do not divide by zero
assert( all(err < absTol) )

figure, hold on, title('power flux density px at boundary')
plot(datRef.w(:)/2/pi, real(pRef(:)), 'x')
plot(dat.w(:)/2/pi, real(p(:)), '.')

%% kinetik energy 
Ek = energyKinetic(gew, dat);
err = (abs(Ek(:) - EkRef(:)))./EkRef(:);
assert( all(err < relTol) )

figure, hold on, title('total kinetik energy')
plot(datRef.w(:)/2/pi, real(EkRef(:)), 'x')
plot(dat.w(:)/2/pi, real(Ek(:)), '.')

%% elastic energy:
Es = energyElastic(gew, dat);
err = (abs(Es(:) - EsRef(:)))./EsRef(:);
assert( all(err < relTol) )

figure, hold on, title('total elastic energy')
plot(datRef.w(:)/2/pi, abs(EsRef(:)), 'x');
plot(dat.w(:)/2/pi, abs(Es(:)), '.');

%% power flux total
P = powerFlux(gew, dat);
err = (abs(P(:) - PRef(:)))./PRef(:);
assert( all(err < relTol) )

figure, hold on, title('total power flux P')
plot(datRef.w(:)/2/pi, real(PRef(:)), 'x')
plot(dat.w(:)/2/pi, real(P(:)), '.')

%% energy vel
ce = energyVel(gew, dat);
err = (abs(ce(:) - ceRef(:)))./ceRef(:);
assert( all(err < relTol) )

figure, hold on, title('energy velocity ce')
plot(datRef.w(:)/2/pi, ceRef(:), 'x')
plot(dat.w(:)/2/pi, ce(:), '.')
