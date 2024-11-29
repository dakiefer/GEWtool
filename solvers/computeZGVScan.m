function [dat] = computeZGVScan(gew, wmax, opts)
% computeZGVScan - Compute ZGV points via iterative shift and search.
%
% Determines zero-group-velocity (ZGV) points (k, w) on the dispersion curves.
% The method rather reliably finds all ZGV points without compromising in speed.
% For this end it first computes candidate wavenumbers close to a shift k0 and
% then refines these with ZGVNewtonBeta(). By iteratively updating the shift k0,
% a (normalized) wavenumber range [kmin, kEnd] is searched for the ZGV points.
% The method is implemented in ZGV_Sylv_MFRDScan() and ZGV_MFRDScan().
% 
% Usage: 
% dat = computeZGVScan(gew):   Returns the ZGV points given the
% guided-wave-description object "gew". ZGV points are searched only up to the 
% 10th cutoff frequency. It is desirable to pass wmax as a second argument
% explicitly, as this might speed up the computation. 
%
% dat = computeZGVScan(gew, wmax):   Also specify the maximum angular frequency
% wmax in rad/s that you are interested in. wmax defines kMax and nModes via the
% longitudinal wave speed. We conjecture that ZGVs cannot exist at phase
% velocities below the minimum bulk longitudinal wave speed.
%
% dat = computeZGVScan(gew, wmax, opts):  Additinally control the algorithm behavior
% with the structure "opts". If you want to use the default value of wmax, pass
% "inf" as its value. The optsions you can specify are (e.g., opts.kEnd=20e3):
% 
% - options in opts:
%      - kEnd: Maximum wavenumber in rad/m to scan to. Reduce kEnd to speed up
%        the computation. ZGV points with kZGV >~ kEnd will not be found. Always 
%        specify kEnd if possible, as the speedup can be considerable.
%      - Neigs: Number of eigenvalues to compute using eigs for a given shift.
%        Reducing Neigs can considerably speed up the computation, especially if
%        you have large matrices. If it is too small, you will miss solutions.
%        The default value is determined depending on the number of modes within the
%        frequency range of interest (wmax). 
%      - Dk: (kMax/(6*nModes+1)): Increment for the target k0 in rad/m where to search
%        for ZGV points. Increase Dk to speed up the computation. Smaller Dk make it
%        more likely that all ZGV points will be found. The default is estimated
%        based on kMax (as for Neigs) and the number of modes in the frequency
%        range of interest (wmax).
%      - kStart: (2*Dk): Initial shift for wavenumber k in rad/m. Attention: searching
%        close to 0 is slow. Don't set it too small.
%      - DeltaPert: (1e-3): Regularization parameter. At ZGV points, we have two
%        very close solutions k and (1+DeltaPert)*k. These are computed as
%        initial guess where the ZGV point is finally computed. Setting
%        DeltaPert too small will lead to an error being thrown by sylvester().
%      - MaxIter: (100) Break search when reaching this number of iterations.
%      - show: (false) whether to display information during calculation 
% 
% For details refer to:
% D. A. Kiefer, B. Plestenjak, H. Gravenkamp, and C. Prada, "Computing 
% zero-group-velocity points in anisotropic elastic waveguides: globally and locally 
% convergent methods." arXiv, Nov. 2022. doi: 10.48550/arXiv.2211.01995.
%
% See also computeZGV, computeZGVDirect, ZGV_Sylv_MFRDScan, ZGV_MFRDScan, ZGVNewtonBeta, Waveguide.
% 
% Design of ZGV_Slyv_MFRDScan() and ZGV_MFRDScan(): B. Plestenjak, University of Ljubljana, Slovenia
% 2022-2023 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France

if nargin < 2, wmax = inf; end
if nargin < 3, opts = struct(); end

% more precise and more beautiful but only since MATLAB R2019b:
% arguments
%     gew           Waveguide
%     wmax    (1,1) double = inf;
%     opts    (1,1) struct = struct();
% end

if ~isscalar(gew)
    compute = @(gewObj) computeZGVScan(gewObj, wmax, opts);
    dat = arrayfun(compute,gew); % apply to every object in the array "gew"
    return; 
end

if isinf(wmax)                  % determine a default maximum frequency 
    wCutoff = gew.cutoffFreq(); % calculate cutoff frequencies
    if length(wCutoff) >= 10
        wmax = wCutoff(10);
    else
        wmax = wCutoff(end);
    end
    warning('GEWTOOL:computeZGV:wmax', 'You did not provide the maximum angular frequency "wmax" you are interested in. I chose %g Hz for you.', wmax/2/pi);
end
% restrict wavenumber search domain according to wmax:
for l = 1:length(gew.lay), mats(l) = gew.lay{l}.mat; end
cList = cell2mat(arrayfun(@(x) x.wavespeeds(), mats, 'UniformOutput', false)); % wave speeds in x-direction
clmin = min(cList(1,:)); % We conjecture that ZGVs cannot exist at phase velocities below the minimum bulk longitudinal wave speed.
kMax = wmax/clmin;  % maximum relevant wavenumber corresponding to wmax: no point in searching above kMax
nModes = gew.nModes(wmax); % number of modes below max. frequency: gives an idea of the number of ZGVs that might occur

% algorithm options have been fine-tuned empirically:
% if ~isfield(opts, 'kStart'),   opts.kStart = 1e-2;                 end % initial shift, i.e., search will be done for k > kStart 
if(~isfield(opts, 'kEnd') || kMax <= opts.kEnd), opts.kEnd = kMax; end % there is no point in searching above kMax   
if ~isfield(opts, 'Neigs'),    opts.Neigs = ceil((nModes)/2)+1;  end % number of candidates to compute at each shift k0
if ~isfield(opts, 'Dk'),       opts.Dk = kMax/(6*nModes+1);         end 
if ~isfield(opts,'show'),      opts.show = false;                  end % do not display iteration results
if ~isfield(opts,'wmax'),      opts.wmax = wmax;                  end 
if opts.show, disp(opts); end
% normalize parameters before passing to ZGV_Sylv_MFRDScan:
opts.kEnd = opts.kEnd*gew.np.h0;
if isfield(opts, 'kStart'), opts.kStart = opts.kStart*gew.np.h0; end
opts.Dk = opts.Dk*gew.np.h0; 
opts.wmax = opts.wmax/gew.np.fh0*gew.np.h0; 

L2 = gew.op.L2; L1 = gew.op.L1; L0 = gew.op.L0; M = gew.op.M;

warnStat = warning('query', 'MATLAB:nearlySingularMatrix');
if ~opts.show, warning('off', 'MATLAB:nearlySingularMatrix'); end
[k, w] = ZGV_Sylv_MFRDScan(L2, L1, L0, M, opts);
if strcmp(warnStat.state, 'on'), warning('on', 'MATLAB:nearlySingularMatrix'); end

kzgv = k/gew.np.h0; 
wzgv = w*gew.np.fh0/gew.np.h0; 
ind = wzgv <= wmax;
dat = GEWdat(gew, kzgv(ind), wzgv(ind));

end
