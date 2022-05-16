%% Dispersion calculation in a generally anisotropic elastic plate
% Depends on the DMSUITE toolbox by Weideman and Reddy: 
% https://mathworks.com/matlabcentral/fileexchange/29-dmsuite
% 
% 2022 - Daniel A. Kiefer
% Institut Langevin, Paris, France
% 
% see also:
% [1] D. A. Kiefer, M. Ponschab, S. J. Rupitsch, and M. Mayle, “Calculating the full 
% leaky Lamb wave spectrum with exact fluid interaction,” The Journal of the Acoustical 
% Society of America, vol. 145, no. 6, pp. 3341–3350, Jun. 2019, doi: 10.1121/1.5109399.
%

h = 1e-3;   % thickness in m
N = 15;     % discretization: polynomial order of interpolants
rho = 7900; lbd = 1.1538e11; mu = 7.6923e10; % steel material

% define and normalize parameters:
II = eye(3).*shiftdim(eye(3), -2); % 4th order "unit tensor"
c = lbd*II + mu*(permute(II, [1 3 4 2]) + permute(II, [1 3 2 4])); % stiffness tensor
c0 = c(1,2,1,2); h0 = h; % normalization parameters
rho0 = rho; fh0 = sqrt(c0/rho0); % normalization parameters
rhon = rho/rho0; cn = c/c0;

% relevant material matrices: 
udof = 3; % Lamb and/or SH
cxx = squeeze(cn(1,udof,udof,1));
cxy = squeeze(cn(1,udof,udof,2)); 
cyx = squeeze(cn(2,udof,udof,1));
cyy = squeeze(cn(2,udof,udof,2));
I = eye(size(cxx)); 


%% mesh: 
dom = [0 1];
% [y, w] = legpts(N); %lobpts
[y, w] = chebpts(N);
% D1 = diffmat(N,1,[0 1],'chebkind2'); D2 = diffmat(N,2,[0 1],'chebkind2');
Psi = chebpoly(0:N, dom);
Psid = diff(Psi);
% mesh = Geometry(dom, N, length(udof));
% mesh.y{1} = y;
me = elemM(Psi);
l2 = elemL2(Psi);
l1 = elemL1(Psi, Psid);
l0 = elemL0(Psid);
M = rhon*me;
L2 = cxx*l2;
L1 = (cxy+cyx)*l1;
L0 = cyy*l0;


%% discretization 
% % [~, Dy_dash] = chebdif(N, 2); 
% % D1 = -2*Dy_dash(:,:,1); % differentiation on unit domain
% % D2 = 4*Dy_dash(:,:,2);  % second order derivative
% D1 = diffmat(N,1,[0 1],'chebkind2'); D2 = diffmat(N,2,[0 1],'chebkind2');
% Id = eye(size(D1));     % identity matrix for discretization
% 
% % % define wave operators:
% L2 = kron(cxx, Id); L1 = kron(cxy + cyx, D1); L0 = kron(cyy, D2); 
% M = kron(rhon*I, Id);
% B1 = kron(cyx, Id([1, N], :)); B0 = kron(cyy, D1([1, N], :));
% 
% % incorporate BCs:
% dofBC = [(0:length(udof)-1)*N+1; (1:length(udof))*N]; % [1, N, N+1, 2*N, 2*N+1, 3*N];
% L2(dofBC, :) = 0; L1(dofBC, :) = B1; L0(dofBC, :) = B0; M(dofBC, :) = 0;
% L2(dofBC, :) = []; L1(dofBC, :) = []; L0(dofBC, :) = []; M(dofBC, :) = [];
% L2(:, dofBC) = []; L1(:, dofBC) = []; L0(:, dofBC) = []; M(:, dofBC) = [];

% % solve for frequency and plot:
kh = linspace(1e-2, 15, 300); % wavenumber*thickness 
whn = nan(size(M, 2), length(kh)); tic 
for ii = 1:length(kh)
    kh0 = kh(ii);
    [wh2] = polyeig((1i*kh0)^2*L2 + (1i*kh0)*L1 + L0, M); 
    whn(:,ii) = sqrt(wh2);
end
fh = (whn/2/pi*fh0); fh(fh == 0) = nan;
chron = toc; fprintf('nF: %d, nK: %d, elapsed time: %g, time per point: %g. ms\n', size(fh, 1), size(fh, 2), chron, chron/length(fh(:))*1e3);

% plot wavenumbers:
kkh = kh.*ones(size(fh));
figure(1), hold on, plot(kkh(:)/h/1e3, fh(:)/h/1e6, '.');
xlim([0, 12]), ylim([0, 6]),
xlabel('k in rad/mm'), ylabel('f in MHz'),




%% element mass and stiffness:
function me = elemM(Psi) 
    me = zeros(length(Psi));
    for i = 1:length(Psi)
        for j = 1:length(Psi)
            me(i,j) = sum(Psi(:,i)*Psi(:,j));
        end
    end
end

function le2 = elemL2(Psi) 
    le2 = zeros(length(Psi));
    for i = 1:length(Psi)
        for j = 1:length(Psi)
            le2(i,j) = sum(Psi(:,i)*Psi(:,j));
        end
    end
end

function le1 = elemL1(Psi, Psid) 
    le1 = zeros(size(Psi,2));
    for i = 1:size(Psi,2)
        for j = 1:size(Psi,2)
            le1(i,j) = sum(Psi(:,i)*Psid(:,j) - Psid(:,i)*Psi(:,j));
        end
    end
end

function le0 = elemL0(Psid) 
    le0 = zeros(size(Psid,2));
    for i = 1:size(Psid,2)
        for j = 1:size(Psid,2)
            le0(i,j) = -sum(Psid(:,i)*Psid(:,j));
        end
    end
end

