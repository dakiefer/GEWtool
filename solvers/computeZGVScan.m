function [dat] = computeZGVScan(gew, wmax, opts)
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

if nargin < 3, opts = []; end
if nargin < 2, wmax = inf; end
L2 = gew.op.L2; L1 = gew.op.L1; L0 = gew.op.L0; M = gew.op.M;

% use sparce matrices when they are big
Ndof = size(M,1);
if Ndof>20
    L0 = sparse(L0.*(abs(L0)>1e-12));
    L1 = sparse(L1.*(abs(L1)>1e-12));
    L2 = sparse(L2.*(abs(L2)>1e-12));
    M = sparse(M.*(abs(M)>1e-12));
end

% if maximum frequency is specified, restrict wavenumber search domain:
waveSpeeds = vertcat(gew.lay.mat.wavespeeds); 
cmin = min(waveSpeeds);
kMax = wmax/cmin*gew.np.h0;

% algorithm options have been fine-tuned empirically:
if isfield(opts, 'kStart'),     opts.kStart = opts.kStart*gew.np.h0; end % normalize
if isfield(opts, 'kMax'),       opts.kMax = opts.kMax*gew.np.h0;     end % normalize
if ~isfield(opts, 'MaxPoints'),     opts.MaxPoints = 50;    end     % number of ZGV points to find
if ~isfield(opts, 'MaxIter'),       opts.MaxIter = 20;      end     % break after reaching MaxIter
if ~isfield(opts, 'kStart'),        opts.kStart = 1;        end     % initial shift kStart, search will be done for k > kStart 
if ~isfield(opts, 'kMax'),          opts.kMax = kMax;       end     % maximum wavenumber to stop scanning
if ~isfield(opts, 'ShiftFactor'),   opts.ShiftFactor = 1.1; end     % relative increment for k0
if ~isfield(opts, 'DeltaPert'),     opts.DeltaPert = 1e-6;  end     % regularization parameter 
if ~isfield(opts, 'Neigs'),         opts.Neigs = 8;         end     % number of candidates to compute at each shift k0
if ~isfield(opts, 'show'),          opts.show = false;      end     % display iteration results

warn = warning('query', 'MATLAB:nearlySingularMatrix');
if ~opts.show, warning('off', 'MATLAB:nearlySingularMatrix'); end
[k, w] = ZGV_MFRDScan(L2, L1, L0, M, opts);
if strcmp(warn.state, 'on'), warning('on', 'MATLAB:nearlySingularMatrix'); end

kzgv = k/gew.np.h0; 
wzgv = w*gew.np.fh0/gew.np.h0; 
ind = wzgv <= wmax;
dat.k = kzgv(ind); dat.w = wzgv(ind);

end
