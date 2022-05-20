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
Ni = 2*N; % integration points (needs to be sufficient to integrate square basis functions)
[yi, w] = chebpts(Ni, dom, 2); % lobpts, legpts
Dy = diffmat(Ni, dom);

Psi = chebpoly(0:N-1, dom);
P = Psi(yi,:);
Pd = squeeze(sum(Dy.*shiftdim(P, -1), 2)); % differentiated polynomials

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
function me = elemM(P, w) 
    PtimesP = P.*permute(P,[1 3 2]);
    me = squeeze( sum(w.'.*PtimesP,1) );
end

function le1 = elemK1(P, Pd, w) 
    PtimesPd = P.*permute(Pd,[1 3 2]);
    le1 = squeeze( sum(w.'.*PtimesPd,1) );
end

function g0 = elemG0(Pd, w) 
    PdtimesPd = Pd.*permute(Pd,[1 3 2]);
    g0 = squeeze( sum(-w.'.*PdtimesPd,1) );
end


end

