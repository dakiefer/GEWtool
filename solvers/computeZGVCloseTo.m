function [dat] = computeZGVCloseTo(gew, varargin)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

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
    if isempty(find(abs(kzgv/ki-1) < 1e-12, 1)) && isempty(find(abs(wzgv/wi-1) < 1e-12, 1))  % not yet in list
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
