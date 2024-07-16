function datZGV = computeZGV(gew, varargin)
% computeZGV - Computes ZGV point iteratively from initial guess.
%
% Determines zero-group-velocity (ZGV) points (k, w) on the dispersion curves via a
% Netwon-type iteration implemented in "ZGVNewtonBeta.m". If initial values (k0, w0) 
% are not provided, you should provide a corresponding waveguide solution "dat", from 
% where intial values will be guessed from. 
% 
% dat = computeZGV(gew, dat): Provide guided-wave-description object
% "gew" and dispersion curve solution "dat". computeZGV will compute
% the group velocity dispersion curves. For each zero-crossing of cg(w), a ZGV point 
% is searched. Only ZGV points with k > 0 are kept, i.e., no cutoff frequencies. 
% 
% dat = computeZGV(gew, w0, k0): Provide guided-wave-description object
% "gew" and an initial guess w0 (angular frequency), k0 (wavenumber) for the ZGV point.
% w0 and k0 can be vectors of initial guesses. At each point [w0(i), k0(i)] a
% ZGV point is searched.
% 
% dat = computeZGV(..., opts): Provide an additinoal structure "opts" containing
% parameters that control the algorithm behavior. See details in
% ZGVNewtonBeta.m.
% 
% For details refer to:
% D. A. Kiefer, B. Plestenjak, H. Gravenkamp, and C. Prada, "Computing 
% zero-group-velocity points in anisotropic elastic waveguides: globally and locally 
% convergent methods." arXiv, Nov. 2022. doi: 10.48550/arXiv.2211.01995.
%
% See also computeZGVScan, computeZGVDirect, ZGVNewtonBeta, Waveguide.
%
% 2022-2024 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France

if nargin == 2 || (nargin == 3 && isstruct(varargin{1})) % dispersion data provided instead of initial guess
    dat = varargin{1};
    if ~isscalar(gew) % compute recursively for every waveguide problem 
        compute = @(gewObj,datObj) computeZGV(gewObj, datObj, varargin{2:end}); % function to apply
        datZGV = arrayfun(compute,gew,dat); % apply to every object in the arrays "gew" and "dat"
        return; 
    end
    if isfield(dat, 'cg')
        cg = dat.cg;
    else
        cg = groupVelAxial(gew, dat);
    end
    sigChange = diff(sign(real(cg)),1,1); % detect where the sign changes
    sigChange = [zeros(1,size(sigChange,2)); sigChange]; % correct size to match cg | w | k
    w0 = dat.w(find(sigChange));
    k0 = dat.k(find(sigChange));
    if nargin == 3, opts = varargin{2}; else, opts = []; end
elseif nargin == 3 && ~isstruct(varargin{1}) || nargin == 4 % initial guess (w0, k0) has been provided
    if ~isscalar(gew) % compute recursively for every waveguide problem 
        compute = @(gewObj) computeZGV(gewObj, varargin{:}); % function to apply
        datZGV = arrayfun(compute,gew); % apply to every object in the array "gew"
        return; 
    end
    w0 = varargin{1}(:); % column vector
    k0 = varargin{2}(:); % column vector
    if nargin == 4, opts = varargin{3}; else, opts = []; end
else
    error('GEWTOOL:computeZGV:wrongNumberOfArguments',...
        'Wrong number of input arguments.');
end
L2 = gew.op.L2; L1 = gew.op.L1; L0 = gew.op.L0; M = gew.op.M;

% algorithm options:
if ~isfield(opts, 'beta_corr'), opts.beta_corr = true;      end % algorithm options
if ~isfield(opts, 'show'),      opts.show = false;          end
if ~isfield(opts, 'maxsteps'),  opts.maxsteps = 10;         end
if ~isfield(opts, 'kmin'),      opts.kmin = 1e-6*min(k0)*gew.np.h0;   end % below kmin -> interprete as cutoff
if ~isfield(opts, 'hermitian'), opts.hermitian = [];        end
if isempty(opts.hermitian)
    if ishermitian(L2) && ishermitian(1i*L1) && ishermitian(L0) && ishermitian(M)
        opts.hermitian = true;
    else 
        opts.hermitian = false; 
    end
end

% initialize:
kzgv = nan(length(k0),1);
wzgv = nan(length(w0),1);
uzgv = nan([length(w0), size(M,1)]);
warn = warning('query', 'MATLAB:nearlySingularMatrix');
if ~opts.show, warning('off', 'MATLAB:nearlySingularMatrix'); end
for i=1:numel(w0)
    if isnan(w0(i)) || isnan(k0(i)) || isinf(w0(i)) || isinf(k0(i))
        warning('GEWTOOL:computeZGV:ignoringInitialGuess',...
            'Ignoring NaN or Inf initial guess.');
        continue; % ignore nan and inf entries
    end
    w0i = w0(i)*gew.np.h0/gew.np.fh0;
    k0i = k0(i)*gew.np.h0;
    if opts.hermitian
        [ki,wi,u] = ZGVNewtonBeta(L2, L1, L0, M, k0i, w0i, [], opts);
    else 
        [ki,wi,u,~] = ZGVNewtonComplex(L2, L1, L0, M, k0i, w0i, [], [], opts);
    end
    notInList = isempty(find(abs(kzgv/ki-1) < 1e-10, 1)) && isempty(find(abs(wzgv/wi-1) < 1e-10, 1));
    if notInList % add to list of converged solutions
        kzgv(i) = ki;
        wzgv(i) = wi;
        uzgv(i,:) = u;
    end
end
if strcmp(warn.state, 'on'), warning('on', 'MATLAB:nearlySingularMatrix'); end

if nargin == 2 % remove nan entries when initial guess vector was not provided manually
    ind = ~isnan(wzgv) & ~isnan(kzgv); 
    wzgv = wzgv(ind); kzgv = kzgv(ind);
    [wzgv, ind] = sort(wzgv); 
    kzgv = kzgv(ind);   % sort in frequency
    uzgv = uzgv(ind,:); % sort in frequency
end

% return as structure:
datZGV.k = kzgv/gew.np.h0; 
datZGV.w = wzgv*gew.np.fh0/gew.np.h0; 
datZGV.u = uzgv;

end
