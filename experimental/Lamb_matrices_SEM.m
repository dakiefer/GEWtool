function [guw] = Lamb_matrices_SEM(mat, h, N)

% parameters: 
c = mat.c; rho = mat.rho;
c0 = c(1,2,1,2); h0 = h; % normalization parameters
rho0 = rho; fh0 = sqrt(c0/rho0); % normalization parameters
rhon = rho/rho0; cn = c/c0;

% relevant material matrices: 
udof = 1:2; % Lamb and/or SH
cxx = squeeze(cn(1,udof,udof,1));
cxy = squeeze(cn(1,udof,udof,2));
cyx = squeeze(cn(2,udof,udof,1));
cyy = squeeze(cn(2,udof,udof,2));
I = eye(size(cxx)); 

%% discretize: 
dom = [0 1];
[y, w] = chebpts(2*N-1, dom, 2); % lobpts, legpts
% D1 = diffmat(N,1,[0 1],'chebkind2'); D2 = diffmat(N,2,[0 1],'chebkind2');
% mesh = Geometry(dom, N, length(udof));
% mesh.y{1} = y;

Psi = chebpoly(0:N-1, dom);
P = zeros(length(y), size(Psi,2));
for i = 1:size(Psi, 2)
    Pi = Psi(:,i);
    P(:,i) = Pi(y);
end

Psid = diff(Psi);
Pd = zeros(length(y), size(Psid,2));
for i = 1:size(Psid, 2)
    Pdi = Psid(:,i);
    Pd(:,i) = Pdi(y);
end

me = elemM(P,w);
k2 = me;
k1 = elemK1(P, Pd, w);
g1 = -k1.';
g0 = elemG0(Pd, w);
M  = kron(rhon*I,me);
K2 = kron(cxx, k2);
K1 = kron(cxy, k1);
G1 = kron(cyx, g1);
G0 = kron(cyy, g0);

L2 = K2; L1 = K1 + G1; L0 = G0;

guw.op.L0 = L0; 
guw.op.L1 = L1;
guw.op.L2 = L2;
guw.op.M = M;
guw.np.fh0 = fh0;
guw.np.h0 = h0;
guw.np.c0 = c0;
guw.np.rho0 = rho0;
guw.geom = Geometry([0, h],N,2);


%% element matrices:
function me = elemM(Psi, w) 
    me = zeros(size(Psi,2));
    for i = 1:size(Psi,2)
        for j = i:size(Psi,2)
            me(i,j) = w*(Psi(:,i).*Psi(:,j));
            me(j,i) = me(i,j);
        end
    end
end

function le1 = elemK1(Psi, Psid, w) 
    le1 = zeros(size(Psi,2));
    for i = 1:size(Psi,2)
        for j = 1:size(Psi,2)
            le1(i,j) = w*(Psi(:,i).*Psid(:,j)); 
        end
    end
end

function g0 = elemG0(Psid, w) 
    g0 = zeros(size(Psid,2));
    for i = 1:size(Psid,2)
        for j = i:size(Psid,2)
            g0(i,j) = -w*(Psid(:,i).*Psid(:,j));
            g0(j,i) = g0(i,j);
        end
    end
end


end

