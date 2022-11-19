% Test orthogonality relation for modes kn at fixed frequency w.
% 
% 2022, Daniel Kiefer, Institut Langevin, ESPCI Paris
%

%% calculate complex wavenumber spectrum at given frequency:
mat = Material('aluminum');
h = 1;
N = 30; 
fh = 8e3;
kh = 1;

plate = Plate(mat, h, N);
gew = plate.Lamb;
dat = computeK(gew,2*pi*fh/h);

%% Test the Auld orthogonality relation for modes kn at fixed frequency w.
% 
% [1] B. A. Auld, Acoustic Fields and Waves in Solids 2, 2nd ed., vol. 2, 2 vols. 
% Malabar: Krieger Publishing Company, 1990.

ind = real(dat.k) >= 0 & imag(dat.k) >= -inf; 
% P = powerFlux(gew, dat); u = dat.u{1}; dat.u{1} = u./sqrt(P); % normalize to unit P
v = velocity(dat); v = v{1}; v = v(1,ind,:,:);
T = stress(gew, dat); T = T{1}; 
tx = permute(T(:,ind,:,1,:), [1 2 3 5 4]);

vm = permute(v, [2 1 3 4]); vn = v;
txm = permute(tx, [2 1 3 4]); txn = tx;

I = sum(-conj(vn).*txm - vm.*conj(txn), 4);
Pmn = 1/4*GEWintegrate(gew, {I}, 3);
Pmn = Pmn/max(abs(Pmn(:))); % normalize

heatmap(log10(abs(Pmn)))
caxis([-8 0]), %xlim([1, 2*N]), ylim([1, 2*N])
title('Auld and Kino orthogonality relation')
% NOTES: 
% - purely imaginary k: single Jordan block
% - purely real k:      diagonal entries
% - complex k:          double Jordan block


%% Test the Fraser orthogonality relation for modes kn at fixed frequency w.
% 
% [1] W. B. Fraser, "Orthogonality relation for the Rayleigh–Lamb modes of vibration 
% of a plate," The Journal of the Acoustical Society of America, vol. 59, no. 1, 
% pp. 215–216, Jan. 1976, doi: 10.1121/1.380851.
% [2] D. Royer and T. Valier-Brasier, Ondes élastiques dans les solides 1: propagation 
% (Elastic waves in solids 1: propagation), vol. 1, 2 vols. London: ISTE éditions, 2021.
% 

ind = real(dat.k) >= 0 & imag(dat.k) >= 0; % does not working for complex pairs?
u = dat.u{1}; u = u(1,ind,:,:);
T = stress(gew, dat); T = T{1}; 
tx = permute(T(:,ind,:,1,:), [1 2 3 5 4]);

uxm = permute(u(:,:,:,1), [2 1 3 4]); uyn = u(:,:,:,2);
txym = permute(tx(:,:,:,2), [2 1 3 4]); txxn = tx(:,:,:,1);

I = uyn.*txym - uxm.*txxn;
Dmn = GEWintegrate(gew, {I}, 3);
Dmn = Dmn(1:end-4, 1:end-4);
Dmnabs = abs(Dmn)/max(abs(Dmn(:)));

heatmap(log10(Dmnabs))
caxis([-8 0]), %xlim([1, 2*N]), ylim([1, 2*N])
title('Fraser orthogonality relation')
% NOTES:
% - with N=30 and fh=1e3 it is accurate for 11 modes -> better accuracy with SEM
% - Seems to assume that k is purely real and bigger than zero


