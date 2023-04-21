function [dat] = computeZGVScan(gew, opts)
% computeZGVScan - Compute ZGV points via iterative shift and search.
%
% Determines zero-group-velocity (ZGV) points (k, w) on the dispersion curves.
% The method rather reliably finds all ZGV points without compromising in speed.
% For this end it first computes candidate wavenumbers close to a shift k0 and
% then refines these with ZGVNewtonBeta(). By iteratively updating the shift k0,
% a (normalized) wavenumber range [1e-8, kmax] is searched for the ZGV points.
% The method is implemented in ZGV_MFRDScan().
% 
% Usage: 
% dat = computeZGVScan(gew):   Returns the ZGV points given the
% guided-wave-description object "gew". 
%
% dat = computeZGVScan(gew, opts):  Additinally control the algorithm behavior
% with the structure "opts". The options are described in ZGV_MFRDScan.m.
% 
% For details refer to:
% D. A. Kiefer, B. Plestenjak, H. Gravenkamp, and C. Prada, "Computing 
% zero-group-velocity points in anisotropic elastic waveguides: globally and locally 
% convergent methods." arXiv, Nov. 2022. doi: 10.48550/arXiv.2211.01995.
%
% See also computeZGV, computeZGVDirect, ZGV_MFRDScan, ZGVNewtonBeta, Waveguide.
% 
% Design of ZGV_MFRDScan(): B. Plestenjak, University of Ljubljana, Slovenia
% 2022 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France

if nargin < 2, opts = []; end
L2 = gew.op.L2; L1 = gew.op.L1; L0 = gew.op.L0; M = gew.op.M;

% % use sparce matrices: faster and more stable because 
% % Rayleigh quotient for mu has a larger imaginary part for full matrices
L0 = sparse(L0);
L1 = sparse(L1);
L2 = sparse(L2);
M = sparse(M);

% algorithm options have been fine-tuned empirically:
if isfield(opts, 'kStart'),     opts.kStart = opts.kStart*gew.np.h0; end % normalize
if isfield(opts, 'kMax'),       opts.kMax = opts.kMax*gew.np.h0;     end % normalize
if ~isfield(opts, 'MaxPoints'),     opts.MaxPoints = 50;    end     % number of ZGV points to find
if ~isfield(opts, 'ShiftFactor'),   opts.ShiftFactor = 1.1; end     % relative increment for k0
if ~isfield(opts, 'DeltaPert'),     opts.DeltaPert = 1e-6;  end     % regularization parameter 
if ~isfield(opts, 'Neigs'),         opts.Neigs = 8;         end     % number of candidates to compute at each shift k0
if ~isfield(opts, 'show'),          opts.show = false;      end     % display iteration results
if ~isfield(opts, 'wmax'),          opts.wmax = inf;        end     % maximum frequency that defines kMax
if ~isfield(opts, 'kStart'),        opts.kStart = 1;        end     % initial shift kStart, search will be done for k > kStart 
% if maximum frequency is specified, restrict wavenumber search domain:
lays = [gew.lay]; mats = [lays.mat];
cList = cell2mat(arrayfun(@(x) x.wavespeeds(), mats, 'UniformOutput', false)); % wave speeds in x-direction
clmin = min(cList(1,:)); % We conjecture that ZGVs cannot exist at phase velocities below the minimum bulk longitudinal wave speed.
kMax = opts.wmax/clmin*gew.np.h0; % to be used only if not provided in opts.kMax
if ~isfield(opts, 'kMax'), opts.kMax = kMax;       end     % maximum wavenumber to stop scanning
if ~isfield(opts, 'MaxIter')    % break after reaching MaxIter
    if isinf(opts.kMax)
        opts.MaxIter = 30;
    else
        opts.MaxIter = inf;
    end
end     

warn = warning('query', 'MATLAB:nearlySingularMatrix');
if ~opts.show, warning('off', 'MATLAB:nearlySingularMatrix'); end
[k, w] = ZGV_MFRDScan(L2, L1, L0, M, opts);
if strcmp(warn.state, 'on'), warning('on', 'MATLAB:nearlySingularMatrix'); end

kzgv = k/gew.np.h0; 
wzgv = w*gew.np.fh0/gew.np.h0; 
ind = wzgv <= opts.wmax;
dat.k = kzgv(ind); dat.w = wzgv(ind);

end
