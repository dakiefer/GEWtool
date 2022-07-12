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
for i=1:numel(w0)
    w0i = w0(i)*gew.np.h0/gew.np.fh0;
    k0i = k0(i)*gew.np.h0;
    mu = w0i^2;
    lambda = 1i*k0i;
    [lambdaZGV,muZGV,u] = ZGV_Lamb_Newton_Hermitian(L0, L1, L2, M, lambda, mu);
    kzgv(i) = -1i*lambdaZGV/gew.np.h0;
    wzgv(i) = sqrt(muZGV)*gew.np.fh0/gew.np.h0;
    [row, col] = ind2sub(size(wzgv), i);
    uzgv(row, col, :) = u;
end

dat.k = real(kzgv); dat.w = real(wzgv); dat.u = uzgv;

end
