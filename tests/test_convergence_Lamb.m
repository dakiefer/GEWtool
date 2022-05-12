%% plot the convergence with respect to the increase in degrees of freedom (dof)
% The convergence is tested against the Rayleigh-Lamb root, i.e., valid for Lamb waves
% in an isotropic plate. The provided approximate solution (w0, k0) is for an
% aluminum plate.

mat = Material('aluminum');
h = 1; 
w0 = 2*pi*5000; % for reference and mode selection
k0 = 5.59596;  % for reference and mode selection (also: 3.6180, 3.8674, 5.59596, 7.8648)
N = 4:30;

%% computing the frequency w
ws=nan(numel(N),1); % allocate
dofs=ws; % allocate
for i=1:numel(N)
    plate = Plate(mat, h, N(i));
    guw = plate.Lamb;
    dofs(i)=size(guw.op.L0,1);
    dat = computeW(guw, k0);
    [~, indSel] = min(abs(dat.w - w0));
    ws(i)=dat.w(indSel);
end

% % compute Rayleigh-Lamb root:
rayLamb = @(wh) abs(rayleighLambS(mat, wh, k0).*rayleighLambA(mat, wh, k0));
options = optimset('Display','iter');
options.TolX = eps;
detRoot = fminbnd(rayLamb, w0-1, w0+1, options);

err=abs(ws-detRoot)/detRoot;

figure, plot(dofs,abs(err),'o--');
ax=gca; ax.YScale='log';
xlabel('matrix size (= 2N)'), ylabel('rel. error')
title('Error w.r.t Rayleigh-Lamb root')

% % test against Rayleigh-Lamb equations:
% errWRef = abs(w0/detRoot-1)
errW = abs(ws(end)/detRoot-1)

%% computing the wavenumbers k
ks=nan(numel(N),1); % allocate
dofs=ks; % allocate
for i=1:numel(N)
    plate = Plate(mat, h, N(i));
    guw = plate.Lamb;
    dofs(i)=size(guw.op.L0,1);
    dat = computeK(guw, w0);
    [~, indSel] = min(abs(dat.k - k0));
    ks(i)=dat.k(indSel);
end

% % compute Rayleigh-Lamb root:
rayLamb = @(kh) abs(rayleighLambS(mat, w0, kh).*rayleighLambA(mat, w0, kh));
options = optimset('Display','iter');
options.TolX = eps;
detRoot = fminbnd(rayLamb, k0-1e-3, k0+1e-3, options);

err=abs(ks-detRoot)/detRoot;

hold on, plot(dofs,abs(err),'o--');
legend({'compute freq w', 'compute wavenumber k'})

% % test against Rayleigh-Lamb equations:
% errKRef = abs(k0/detRoot-1)
errK = abs(ks(end)/detRoot-1)

