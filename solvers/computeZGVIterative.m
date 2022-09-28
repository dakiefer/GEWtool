function [dat] = computeZGVIterative(gew, varargin)
% computeZGVIterative - Computes ZGV points.
% Determines zero-group-velocity (ZGV) points (k, w) on the dispersion curves via a
% Netwon-type iteration implemented in "ZGV_NewtonBeta.m". If initial values (k0, w0) 
% are not provided, you should provide a corresponding waveguide solution "dat", from 
% where intial values will be guessed from. 
% 
% dat = computeZGVIterative(gew, dat): Provide waveguide description object
% "gew" and dispersion curve solution "dat". computeZGVIterative will compute
% the group velocity dispersion curves. For each zero-crossing, a ZGV point is
% searched. Only true ZGV points with k > 0 are kept. 
% 
% dat = computeZGVIterative(gew, w0, k0): Provide waveguide description object
% "gew" and an initial guess w0 (angular frequency), k0 (wavenumber) for the ZGV point.
% w0 and k0 can be vectors of initial guesses. At each point [w0(i), k0(i)] a
% ZGV point is searched.
%
% See also computeZGVDirect, ZGV_NewtonBeta, Waveguide.
% 
% 2022 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France
% 

if nargin == 2 % initial guess (w0, k0) where cg changes sign
    dat = varargin{1};
    cg = groupVel(gew, dat);
    sigChange = diff(sign(real(cg)),1,2); % detect where the sign changes
    w0 = dat.w(find(sigChange));
    k0 = dat.k(find(sigChange));
elseif nargin == 3 % initial guess (w0, k0) has been provided
    w0 = varargin{1};
    k0 = varargin{2};
else
    error('GEWTOOL:computeZGVCloseTo:wrongNumberOfArguments', 'Wrong number of input arguments.');
end

L2 = gew.op.L2; L1 = gew.op.L1; L0 = gew.op.L0; M = gew.op.M;

kzgv = nan(size(k0));
wzgv = nan(size(w0));
uzgv = nan([size(w0), size(M,1)]);
opts.beta_corr = true; % algorithm options
opts.show = false;
opts.maxsteps = 10;
w = warning('query', 'MATLAB:nearlySingularMatrix');
warning('off', 'MATLAB:nearlySingularMatrix')
for i=1:numel(w0)
    if isnan(w0(i)) || isnan(k0(i)) || isinf(w0(i)) || isinf(k0(i))
        warning('GEWTOOL:computeZGVCloseTo:ignoringInitialGuess', 'Ignoring NaN or Inf initial guess.');
        continue; % ignore nan and inf entries
    end
    w0i = w0(i)*gew.np.h0/gew.np.fh0;
    k0i = k0(i)*gew.np.h0;
    [ki,wi,u] = ZGV_NewtonBeta(L0, L1, L2, M, k0i, w0i, [], opts); %wi = sqrt(mui);
%     ki = ki/gew.np.h0; wi = sqrt(mui)*gew.np.fh0/gew.np.h0;
    isCutOff = (ki/k0i) < 1e-12;
    notInList = isempty(find(abs(kzgv/ki-1) < 1e-12, 1)) && isempty(find(abs(wzgv/wi-1) < 1e-12, 1));  % not yet in list
    if ~isCutOff && notInList
        kzgv(i) = ki;
        wzgv(i) = wi;
        [row, col] = ind2sub(size(wzgv), i);
        uzgv(row, col, :) = u;
    end
end
if strcmp(w.state, 'on'), warning('on', 'MATLAB:nearlySingularMatrix'); end

dat.k = real(kzgv)/gew.np.h0; 
dat.w = real(wzgv)*gew.np.fh0/gew.np.h0; 
dat.u = uzgv;

end
