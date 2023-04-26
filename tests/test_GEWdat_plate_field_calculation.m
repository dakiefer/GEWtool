% Run using: runtests()
% Test preprocessing implemented in GEWdat/ to pre-computed results.
%
% see also: 
% https://fr.mathworks.com/help/matlab/matlab_prog/write-script-based-unit-tests.html
%
% 2022 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France

% % options : 
show = false;

% % load pre-computed fields:
load('data/Lamb_ref_fields.mat')
warnStat = warning; warning off; % save state and turn off

% % current implementation:
N = 30; % needs to be high enough such that also the evanescent waves converge
plate = Plate(mat, h, N); % mat and h from Lamb_ref_fields.mat
gew = plate.Lamb;
% gew = gew.linearizeInK;
clear opts;
opts.parallel = false;
opts.subspace = false;
opts.sparse = false;
dat = computeK(gew, 2*pi*freq, Nmodes, opts); % freq and Nmodes from Lamb_ref_fields.mat
[y, wInt] = Layer.nodes(N); wInt = wInt*h; y = (y-0.5)*h;
dat.u{1} = normalizeL2(dat.u{1}, wInt); % % normalize to ∫ conj(u).u dy = 1
relTol = 1e-4;

%% displacements:
% No interpolation needed, as we compare at the boundary.
u = abs(dat.u{1}(:,:,1,1));
err = (abs(u(:) - uRef(:)))./uRef(:);
assert( all(err < relTol) )

if exist('show', 'var') && show
    figure, hold on, title('displacement ux at boundary')
    plot(datRef.w(:)/2/pi, abs(uRef(:)), 'x')
    plot(dat.w(:)/2/pi, abs(u(:)), '.')
end

%% particle velocity:
v= velocity(dat);
v = abs(v{1});
v = v(:,:,1,1);
err = (abs(v(:) - vRef(:)))./vRef(:);
assert( all(err < relTol) )

if exist('show', 'var') && show
    figure, hold on, title('velocity vx at boundary')
    plot(datRef.w(:)/2/pi, real(vRef(:)), 'x')
    plot(dat.w(:)/2/pi, real(v(:)), '.')
end

%% strain: 
S = strain(gew, dat);
S = abs(S{1});
S = S(:,:,1,1,1);
err = (abs(S(:) - SRef(:)))./SRef(:);
assert( all(err < relTol) )

if exist('show', 'var') && show
    figure, hold on, title('strain Sxx at boundary')
    plot(datRef.w(:)/2/pi, real(SRef(:)), 'x')
    plot(dat.w(:)/2/pi, S(:), '.')
end

%% stress: 
T0 = 1; % mean(TRef(:));
TRefN = TRef/T0;
T = stress(gew, dat);
T = abs(T{1}(:,:,1,1,1))/T0; % complex-valued, Txx at top is generally non-zero
err = (abs(T(:) - TRefN(:)))./TRefN(:);
assert( all(err < relTol) )

if exist('show', 'var') && show
    figure, hold on, title('stress Txx at boundary')
    plot(datRef.w(:)/2/pi, TRefN(:), 'x')
    plot(dat.w(:)/2/pi, T(:), '.'), %ylim([0, 3e14])
end

%% power flux at interface 
% modes are sorted by wavenumber magnitude, i.e., modes with positive and
% negative real part of the wavenumber are in arbitrary order. To ensure that
% the correct modes are compared to each other, we compare only the ones with
% positive real wavenumber.
indRef = real(datRef.k) > 0;
ind =    real(dat.k) > 0;
pRefPos = pRef(indRef); % only with real(k) > 0
p = poyntingVec(gew, dat);
p = p{1}(:,:,1,1);
p = p(ind); % only with real(k) > 0
err = abs(p - pRefPos)./pRefPos; % relative error
assert( all(err < relTol) )

if exist('show', 'var') && show
    figure, hold on, title('power flux density px at boundary')
    plot(datRef.w(indRef)/2/pi, real(pRefPos), 'x')
    plot(dat.w(ind)/2/pi, real(p), '.')
end

%% kinetik energy 
Ek = energyKinetic(gew, dat);
err = (abs(Ek(:) - EkRef(:)))./EkRef(:);
assert( all(err < relTol) )

if exist('show', 'var') && show
    figure, hold on, title('total kinetik energy')
    plot(datRef.w(:)/2/pi, real(EkRef(:)), 'x')
    plot(dat.w(:)/2/pi, real(Ek(:)), '.')
end

%% elastic energy:
Es = energyElastic(gew, dat);
err = (abs(Es(:) - EsRef(:)))./EsRef(:);
assert( all(err < relTol) )

if exist('show', 'var') && show
    figure, hold on, title('total elastic energy')
    plot(datRef.w(:)/2/pi, abs(EsRef(:)), 'x');
    plot(dat.w(:)/2/pi, abs(Es(:)), '.');
end

%% power flux total
% modes are sorted by wavenumber magnitude, i.e., modes with positive and
% negative real part of the wavenumber are in arbitrary order. To ensure that
% the correct modes are compared to each other, we compare only the ones with
% positive real wavenumber.
indRef = real(datRef.k) > 0;
ind =    real(dat.k) > 0;
PRefPos = PRef(indRef);
P = powerFlux(gew, dat);
PPos = P(ind);
err = (abs(PPos - PRefPos))/mean(abs(PRefPos)); % absolute normalized error
assert( all(err < relTol) )

if exist('show', 'var') && show
    figure, hold on, title('total power flux P')
    plot(datRef.w(indRef)/2/pi, real(PRefPos(:)), 'x')
    plot(dat.w(ind)/2/pi, real(PPos(:)), '.')
end

%% energy vel
% modes are sorted by wavenumber magnitude, i.e., modes with positive and
% negative real part of the wavenumber are in arbitrary order. To ensure that
% the correct modes are compared to each other, we compare only the ones with
% positive real wavenumber.
indRef = real(datRef.k) > 0;
ind =    real(dat.k) > 0;
ceRefPos = ceRef(indRef);
ce = energyVel(gew, dat);
cePos = ce(ind);
err = (abs(cePos - ceRefPos))/mean(abs(ceRefPos)); % absolute normalized error
assert( all(err < relTol) )

if exist('show', 'var') && show
    figure, hold on, title('energy velocity ce')
    plot(datRef.w(indRef)/2/pi, ceRefPos, 'x')
    plot(dat.w(ind)/2/pi, cePos, '.')
end

%% reset state of warnings
% new block is run even if the previous tests fail (using "runtests")
warning(warnStat); % reset warning state
