function [k,w,u,flag,err] = ZGV_NewtonBeta(L0, L1, L2, M, k, w, u, opts)

% [k,w,u,flag,err] = ZGV_NewtonBeta(L0,L1,L2,M,k,mu,u,opts) 
% Returns ZGV points for a twoparameter pencil 
% (L0 + 1i*k*L1 + 1i*k^2*L2 + w^2*M)u = 0 
% using the Newton method with complex correction. It is assumed that 
% A = L0 + 1i*k*L1 - k^2*L2 + w^2*M is hermitian. This is the case
% for a Galerkin discretization (e.g. FE) of a lossless waveguide.
% 
% Output: 
%    - a point (k,w) and unit vector u such that:
%      a) (L0 + 1i*k*L1 - k^2*L2 + w^2*M)*u = 0  
%      b) u'*(1i*L1 - 2*k*L2)*u = 0
%    - flag: convergence (true) or not (false)
%    - err: norm of the residual
%
% Input: 
%    - L0, L1, L2, M: square matrices such that L0, 1i*L1 and L2 are hermitian
%    - k, w: initial approximation for the ZGV point
%    - u: initial approximation for the right eigenvector, if
%      not provided or empty then the singular vectors corresponding to the 
%      smallest singular value of (L0 + 1i*k*L1 + 1i*k^2*L2 + w^2*M) are used
%    - options in opts:
%         - maxsteps: (10) maximal number of steps 
%         - tol: (1e-14) tolerance (relative to the norm of matrices)
%         - show: (false), set to 'true' to display values and residuals
%
% 2022 Bor Plestenjak, Daniel A. Kiefer

narginchk(6, 8);
if nargin < 7, u = [];    end
if nargin < 8, opts = []; end

if isfield(opts,'show'),      show     = opts.show;       else, show     = false;                       end
if isfield(opts,'tol'),       tol      = opts.tol;        else, tol      = 1e-14;                   end
if isfield(opts,'maxsteps'),  maxsteps = opts.maxsteps;   else, maxsteps = 10;                      end
if isfield(opts,'beta_corr'), beta_corr = opts.beta_corr; else, beta_corr = true;                      end

mu = w^2; 
if isempty(u)
    [~,~,V] = svd(L0 + k*1i*L1 - k^2*L2 + mu*M);
    u = V(:,end);
end

% we use relative tolerance
reltol = tol * norm([L0 1i*L1 -L2 M],'fro');

n = size(L0,1);
step = 0;
flag = false;
tmp = zeros(n+2,1); 
tmp(n+1,1) = 1;
while step < maxsteps
    step = step + 1;
    A = L0 + k*1i*L1 - k^2*L2 + mu*M;
    B = 1i*L1 - 2*k*L2;
    J = [A            B*u         M*u;
         2*u'*B     -2*u'*L2*u   0;
         2*u'           0         0   ];
    res = [ A*u; u'*B*u; u'*u-1 ];
    err = norm(res);
    if show
        fprintf('step %3d | err: %8.3e | lambda (%16.9e,%16.9e) | mu (%16.9e,%16.9e) \n',step,err,real(k),imag(k),real(mu),imag(mu));
    end
    if err<reltol
        flag = true;
        break
    end
    if beta_corr
        kor = -J\[res tmp];
        p = kor(:,1);
        q = kor(:,2);
        beta = -imag(p(n+1,1))/real(q(n+1,1)); % correction suggested by Lu & Su
        if abs(beta)>1e-12
%             fprintf('     correction beta : %8.3e\n',beta)
            upd = p + 1i*beta*q;
        else
            upd = p;
        end
    else
        upd = -J\res;
    end
    u = u + upd(1:n);
    k = k + real(upd(n+1));
    mu = mu + real(upd(n+2));
    u = u/norm(u);
end

w = real(sqrt(mu));
if ~flag
    k = nan; 
    w = nan; 
    u = nan(size(u));
end

