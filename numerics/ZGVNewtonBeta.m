function [k,w,u,isConverged,err] = ZGVNewtonBeta(L2, L1, L0, M, k, w, u, opts)
% ZGVNewtonBeta - Newton-type iteration to locate ZGV points. 
% Computes a ZGV point given an initial guess (w0, k0) and optionally u0. 
% A Newton-type method with complex correction is used for this end [1]. If the
% iteration converges, the solutions (k, w, u) satisfies: 
%   W*u   = [(1i*k)^2*L2 + 1i*k*L1 + L0 + w^2*M]*u = 0,     (waveguide problem)
%   u'*iWd*u = u'*[-2*k*L2 + 1i*L1]*u = 0                      (group velocity is zero)
%   u'*u - 1 = 0                                    (unit eigenvector u)
% The method assumes that (1i*k)^2*L2 + 1i*k*L1 + L0 + w^2*M is Hermitian. This is 
% the case for a Galerkin discretization (e.g. FE) of a lossless waveguide.
% 
% Output: 
%    - ZGV point (k,w) and corresponding unit vector u
%    - flag: convergence (true) or not (false)
%    - err: norm of the residual
%
% Input: 
%    - L2, L1, L0, M: square matrices such that L2, 1i*L1 and L1 are Hermitian
%    - k, w: initial approximation for the ZGV point
%    - u (optional): initial approximation for the right eigenvector, if
%      not provided or empty then the singular vectors corresponding to the 
%      smallest singular value of [(1i*k)^2*L2 + 1i*k*L1 + L0 + w^2*M] are used
%    - options in structure opts:
%         - maxsteps:   (10) maximum number of steps 
%         - tol:        (1e-14) tolerance (relative to the norm of matrices)
%         - beta_corr:  (true) turn on/off complex correction
%         - show:       (false) set to 'true' to display values and residuals
%
% Literature:
% [1] D. A. Kiefer, B. Plestenjak, H. Gravenkamp, and C. Prada, “Computing 
% zero-group-velocity points in anisotropic elastic waveguides: globally and locally 
% convergent methods.” arXiv, Nov. 03, 2022. doi: 10.48550/arXiv.2211.01995.
%
% 2022 Bor Plestenjak, Daniel A. Kiefer

% % check arguments:
narginchk(6, 8);
if nargin < 7, u = [];    end
if nargin < 8, opts = []; end
if isfield(opts,'show'),      show     = opts.show;       else, show     = false;  end
if isfield(opts,'tol'),       tol      = opts.tol;        else, tol      = 1e-14;  end
if isfield(opts,'maxsteps'),  maxsteps = opts.maxsteps;   else, maxsteps = 10;     end
if isfield(opts,'beta_corr'), beta_corr = opts.beta_corr; else, beta_corr = true;  end
if ~isreal(k) || ~isreal(w)
    error('ZGV_NewtonBeta: initial guess (k, w) should be real-valued!');
end

% % initialize: 
mu = w^2; 
if isempty(u) % use singular vector of smallest singular value as initial guess for u
    [~,~,V] = svd(-k^2*L2 + 1i*k*L1 + L0 + mu*M);
    u = V(:,end);
end

reltol = tol * norm([L0 1i*L1 -L2 M],'fro'); % relative tolerance to end iteration

n = size(L0,1);
step = 0;            % step number
isConverged = false; % initialize as false
ek = zeros(n+2,1); ek(n+1,1) = 1; % unit vector such that p.ek = k, p = [u, k, w^2]
while step < maxsteps
    % % check first if the solution is converged (satisfies the above system):
    W  = -k^2*L2 + k*1i*L1 + L0 + mu*M;
    iWd = -2*k*L2 +   1i*L1;
    Fp = [ W*u; u'*iWd*u; u'*u-1 ]; % residuum: F(p) should be zero, with p = [u, k, w^2]
    err = norm(Fp);
    if show
        fprintf('step %3d: |res| = %8.3e,  k = %16.9e,  w = %16.9e \n',...
            step,err,real(k),sqrt(mu));
    end
    if err<reltol
        isConverged = true;
        break
    end
    % % Solution is not yet converged. Perform an update for p = [u, k, w^2]:
    step = step + 1;
    J = [ W            iWd*u         M*u;
          2*u'*iWd   -2*u'*L2*u       0;
          2*u'           0            0   ];
    if beta_corr % update with complex correction
        kor = -J\[Fp ek]; % calculate [q; s]
        q = kor(:,1);
        s = kor(:,2);
        beta = -imag(q(n+1))/real(s(n+1)); % correction factor suggested by Lu & Su
        if abs(beta)>1e-12 % only perform complex correction if beta is non-negligible
            if show, fprintf('     correction beta : %8.3e\n',beta); end
            Delta_p = q + 1i*beta*s; % update for p with complex correction
        else
            Delta_p = q;
        end
    else % updated without complex correction
        Delta_p = -J\Fp;
    end
    u = u + Delta_p(1:n);
    k = k + real(Delta_p(n+1));   % imaginary part is only due to numerical precission
    mu = mu + real(Delta_p(n+2)); % imaginary part is only due to numerical precission
    u = u/norm(u);
end

% return values:
w = real(sqrt(mu)); % angular frequency
if ~isConverged % return nan if not converge, i.e., residuum is larger then tol
    k = nan;
    w = nan;
    u = nan(size(u));
end

