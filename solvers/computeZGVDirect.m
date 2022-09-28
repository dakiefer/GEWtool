function [dat] = computeZGVDirect(gew)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

L2 = gew.op.L2; L1 = gew.op.L1; L0 = gew.op.L0; M = gew.op.M;
L1 = 1i*L1;
L2 = -L2;

opts.MaxPoints = 50;
opts.MaxIter = 20;
opts.ZeroShift = 0.5;

% for N>10 we make matrices sparse
Ndof = size(M,1);
if Ndof>30
    L0 = sparse(L0.*(abs(L0)>1e-12));
    L1 = sparse(L1.*(abs(L1)>1e-12));
    L2 = sparse(L2.*(abs(L2)>1e-12));
    M = sparse(M.*(abs(M)>1e-12));
end

ZGVpoints = scan_ZGV_modes(L0,L1,L2,M,gew.np.fh0/gew.np.h0,opts);

dat.k = real(ZGVpoints(:,1))/gew.np.h0; 
dat.w = 2*pi*real(ZGVpoints(:,2)); 
% dat.u = uzgv;

end