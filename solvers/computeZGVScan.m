function [dat] = computeZGVScan(gew)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

L2 = gew.op.L2; L1 = gew.op.L1; L0 = gew.op.L0; M = gew.op.M;

% sparce matrices when they are big
Ndof = size(M,1);
if Ndof>30
    L0 = sparse(L0.*(abs(L0)>1e-12));
    L1 = sparse(L1.*(abs(L1)>1e-12));
    L2 = sparse(L2.*(abs(L2)>1e-12));
    M = sparse(M.*(abs(M)>1e-12));
end

opts.MaxPoints = 50;
opts.MaxIter = 20;
opts.ZeroShift = 1e-3; % start searching close to k = 0
opts.ShiftFactor = 1.1;
opts.Neigs = 8;
[k, w] = ZGV_MFRDScan(L0,L1,L2,M,opts);

dat.k = k/gew.np.h0; 
dat.w = w*gew.np.fh0/gew.np.h0; 
% dat.u = uzgv;

end