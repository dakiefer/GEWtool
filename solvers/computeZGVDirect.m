function [dat] = computeZGVDirect(gew, opts)
% computeZGVDirect - Compute ZGV points by solving the singular three-parameter EVP.
%
% Determines zero-group-velocity (ZGV) points (k, w) on the dispersion curves.
% The method guarantees to find all ZGV points without the need for initial guesses. 
% This comes at rather high computational 
% cost. The method should only be used for problems with small matrix size (~40x40),
% otherwise use computeZGVScan() or computeZGVIterative(). 
% The method is implemented in ZGVDirect() and relies on the "MultiParEig" Toolbox 
% by Bor Plestenjak. Get version 2.6 or later here:
% https://www.mathworks.com/matlabcentral/fileexchange/47844-multipareig
% 
% For details refer to:
% D. A. Kiefer, B. Plestenjak, H. Gravenkamp, and C. Prada, "Computing 
% zero-group-velocity points in anisotropic elastic waveguides: globally and locally 
% convergent methods." arXiv, Nov. 2022. doi: 10.48550/arXiv.2211.01995.
%
% See also computeZGVIterative, computeZGVScan, ZGVDirect, Waveguide.
% 
% Design of ZGVDirect(): B. Plestenjak, University of Ljubljana, Slovenia
% 2022 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France

L2 = gew.op.L2; L1 = gew.op.L1; L0 = gew.op.L0; M = gew.op.M;

if nargin<2, opts=[]; end
if ~isfield(opts,'sc_steps'),  opts.sc_steps=2;   end
if ~isfield(opts,'showrank'),  opts.showrank=1;   end
if ~isfield(opts,'rrqr'),      opts.rrqr=1;       end
if ~isfield(opts,'membtol'),   opts.membtol=1e-4; end

[k, w] = ZGVDirect(L2,L1,L0,M,opts);

dat.k = k/gew.np.h0; 
dat.w = w*gew.np.fh0/gew.np.h0;
% dat.u = uzgv;

end