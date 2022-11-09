function [dat] = computeZGVScan(gew)
% computeZGVScan - Compute ZGV points via iterative shift and search.
%
% Determines zero-group-velocity (ZGV) points (k, w) on the dispersion curves.
% The method rather reliably finds all ZGV points without compromising in speed.
% For this end it first computes candidate wavenumbers close to a shift k0 and
% then refines these with ZGVNewtonBeta(). By iteratively updating the shift k0,
% a (normalized) wavenumber range [1e-8, kmax] is searched for the ZGV points.
% The method is implemented in ZGV_MFRDScan().
% 
% For details refer to:
% D. A. Kiefer, B. Plestenjak, H. Gravenkamp, and C. Prada, "Computing 
% zero-group-velocity points in anisotropic elastic waveguides: globally and locally 
% convergent methods." arXiv, Nov. 2022. doi: 10.48550/arXiv.2211.01995.
%
% See also computeZGVIterative, computeZGVDirect, ZGV_MFRDScan, ZGVNewtonBeta, Waveguide.
% 
% Design of ZGV_MFRDScan(): B. Plestenjak, University of Ljubljana, Slovenia
% 2022 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France

L2 = gew.op.L2; L1 = gew.op.L1; L0 = gew.op.L0; M = gew.op.M;

% use sparce matrices when they are big
Ndof = size(M,1);
if Ndof>20
    L0 = sparse(L0.*(abs(L0)>1e-12));
    L1 = sparse(L1.*(abs(L1)>1e-12));
    L2 = sparse(L2.*(abs(L2)>1e-12));
    M = sparse(M.*(abs(M)>1e-12));
end

% algorithm options have been fine-tuned empirically:
opts.MaxPoints = 50;     % number of ZGV points to find
opts.MaxIter = 20;       % break after reaching MaxIter
opts.kStart = 1;         % initial shift kStart, search will be done for k > kStart 
opts.ShiftFactor = 1.1;  % relative increment for k0
opts.DeltaPert = 1e-6;   % regularization parameter 
opts.Neigs = 8;          % number of candidates to compute at each shift k0
opts.show = false;       % display iteration results
[k, w] = ZGV_MFRDScan(L2, L1, L0, M, opts);

dat.k = k/gew.np.h0; 
dat.w = w*gew.np.fh0/gew.np.h0; 
% dat.u = uzgv;

end