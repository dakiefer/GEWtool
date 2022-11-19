%% plot the convergence with respect to the increase in degrees of freedom (dof)
% The convergence is tested against the Rayleigh-Lamb root, i.e., valid for Lamb waves
% in an isotropic plate. The provided approximate solution (w0, k0) is for an
% aluminum plate.

clear all
addpath('../experimental')
mat = Material('aluminum');
h = 1; 
w0 = 2*pi*5000; % for reference and mode selection
k0 = 5.59596;  % for reference and mode selection (also: 3.6180, 3.8674, 5.59596, 7.8648)
N = 4:30;

%% computing the frequency w
ws=nan(numel(N),1); % allocate
timings = nan(numel(N),1);
dofs=ws; % allocate
disp('Test computeW():'), 
for i=1:numel(N)
%     plate = Plate(mat, h, N(i)); gew = plate.Lamb;
%     gew = Lamb_matrices_rectangularSCM(mat, h, N(i));
    gew = Lamb_matrices_SEM(mat, h, N(i));
    ftest = @() Lamb_matrices_SEM(mat, h, N(i));
    timings(i) = timeit(ftest);
    dofs(i)=size(gew.op.L0,1);
    dat = computeW(gew, k0);
    [~, indSel] = min(abs(dat.w - w0));
    ws(i)=dat.w(indSel);
end
figure(1), hold on, plot(dofs, timings);
title('timing'), xlabel('dofs'), ylabel('time in s')

% % compute Rayleigh-Lamb root:
rayLamb = @(wh) rayleighLambS(mat, wh, k0).*rayleighLambA(mat, wh, k0);
options = optimset('Display','iter');
options.TolX = eps;
detRoot = fzero(rayLamb, w0, options); % only for real arguments
residuumAtRLRoot = abs(rayLamb(detRoot))

err=abs(ws-detRoot)/detRoot;
errWNmax = err(end)

open('data/convergence_SCM.fig')
figure(2), hold on, plot(dofs,abs(err),'o--');
ax=gca; ax.YScale='log';
xlabel('matrix size (= 2N)'), ylabel('rel. error')
title('Error w.r.t Rayleigh-Lamb root')
legend({'SCM', 'new'})



%% computing the wavenumbers k
% ks=nan(numel(N),1); % allocate
% dofs=ks; % allocate
% disp('Test computeK():'), tic
% for i=1:numel(N)
%     plate = Plate(mat, h, N(i));  gew = plate.Lamb;
% %     gew = Lamb_matrices_rectangularSCM(mat, h, N(i));
%     dofs(i)=size(gew.op.L0,1);
%     dat = computeK(gew, w0);
%     [~, indSel] = min(abs(dat.k - k0));
%     ks(i)=dat.k(indSel);
% end
% toc
% 
% % % compute Rayleigh-Lamb root:
% rayLamb = @(kh) rayleighLambS(mat, w0, kh).*rayleighLambA(mat, w0, kh);
% options = optimset('Display','iter');
% options.TolX = eps;
% detRoot = fzero(rayLamb, k0, options); % only for real arguments
% residuumAtRLRoot = abs(rayLamb(detRoot))
% 
% err=abs(ks-detRoot)/detRoot;
% errKNmax = err(end)
% 
% hold on, plot(dofs,abs(err),'o--');
% legend({'compute freq w', 'compute wavenumber k'})


