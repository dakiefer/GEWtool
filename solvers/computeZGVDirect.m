function [dat] = computeZGVDirect(gew, opts)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

L2 = gew.op.L2; L1 = gew.op.L1; L0 = gew.op.L0; M = gew.op.M;

if nargin<2, opts=[]; end
if ~isfield(opts,'sc_steps'),  opts.sc_steps=2;   end
if ~isfield(opts,'showrank'),  opts.showrank=1;   end
if ~isfield(opts,'rrqr'),      opts.rrqr=1;       end
if ~isfield(opts,'membtol'),   opts.membtol=1e-4; end

[k, w] = ZGVDirect(L0,L1,L2,M,opts);

dat.k = k/gew.np.h0; 
dat.w = w*gew.np.fh0/gew.np.h0;
% dat.u = uzgv;

end