function [k,w,x,y,isConverged,err] = ZGVNewtonComplex(L2,L1,L0,M,k,w,x,y,opts)
% ZGVNewtonComplex - Newton iteration to locate ZGV points of non-Hermitian matrices. 
% Computes a ZGV point given an initial guess (k0, w0) for the wavenumber and
% angular frequency and optionally the right and left eigenvectors x0, y0. A
% Newton iteration is used for the computation. If the iteration converges, the
% solutions (k, w, x, u) satisfies:
%   W*x   = [(1i*k)^2*L2 + 1i*k*L1 + L0 + w^2*M]*x = 0,   (waveguide problem)
%   y'*iWd*x = y'*[-2*k*L2 + 1i*L1]*x = 0                 (group velocity is zero)
%   x'*x - 1 = 0                                          (unit right eigenvector x)
%   y'*y - 1 = 0                                          (unit left eigenvector y)
% Contrary to ZGVNewtonBeta(), this method does not assume that (1i*k)^2*L2 +
% 1i*k*L1 + L0 + w^2*M be Hermitian. Non-Hermitian matrices might be obtained
% even for a lossless waveguide for certain numerical discretization methods,
% e.g., the Spectral Collocation Method.
% 
% Output: 
%    - ZGV point (k,w) and corresponding unit right and left eigenvectors x, y.
%    - isConverged: convergence (true) or not (false)
%    - err: norm of the residual
%
% Input: 
%    - L2, L1, L0, M: waveguide operator matrices
%    - k, w: initial approximation for the ZGV point
%    - x, y (optional): initial approximation for the right and left eigenvectors. 
%      If not provided or empty, then the singular vectors corresponding to the 
%      smallest singular value of [(1i*k)^2*L2 + 1i*k*L1 + L0 + w^2*M] are used.
%    - options in structure opts:
%         - maxsteps:   (10) maximum number of steps 
%         - tol:        (1e-14) tolerance (relative to the norm of matrices)
%         - beta_corr:  (true) turn on/off complex correction
%         - show:       (false) set to 'true' to display values and residuals
%         - kmin:       (1e-6) below this value, the result is interpreted as a
%                       cutoff frequency. 
%         - hermitian:  whether waveguide problem is Hermitian. If not provided,
%                       this is determined from the matrices L2, L1, L0, M.
%
% Literature:
% [1] D. A. Kiefer, B. Plestenjak, H. Gravenkamp, and C. Prada, “Computing 
% zero-group-velocity points in anisotropic elastic waveguides: globally and locally 
% convergent methods.” arXiv, Nov. 03, 2022. doi: 10.48550/arXiv.2211.01995.
%
% 2022-2024 Bor Plestenjak, Daniel A. Kiefer

narginchk(6, 9);
if nargin < 7, x = [];    end
if nargin < 8, y = [];    end
if nargin < 9, opts = []; end
if isfield(opts,'show'),      show     = opts.show;       else, show     = false;  end
if isfield(opts,'tol'),       tol      = opts.tol;        else, tol      = 1e-14;  end
if isfield(opts,'maxsteps'),  maxsteps = opts.maxsteps;   else, maxsteps = 10;     end
if isfield(opts, 'kmin'),     kmin = opts.kmin;           else, kmin = 1e-4*k;     end % below kmin -> interprete as cutoff

mu = w^2;
if isempty(x) || isempty(y)
    [U,~,V] = svd(L0 + 1i*k*L1 - k^2*L2 + mu*M);
    if isempty(x), x = V(:,end); end
    if isempty(y), y = U(:,end); end
end

% we use relative tolerance
reltol = tol * norm([L0 L1 L2 M],'fro');

n = size(L0,1);
step = 0;
isConverged = false;
while step < maxsteps
    step = step + 1;
    y = conj(y); 
    A = L0 + 1i*k*L1 - k^2*L2 + mu*M;
    DA = 1i*L1 - 2*k*L2;
    DAx = DA*x;
    DAy = DA.'*y;
    J = [A          zeros(n)   DAx      M*x;
         zeros(n)   A.'        DAy    M.'*y;
         DAy.'      DAx.'    -2*y.'*L2*x  0;
         x'       zeros(1,n) 0     0;
         zeros(1,n) y'         0     0];
    res = [A*x; A.'*y; y.'*DAx; (x'*x-1)/2; (y'*y-1)/2];
    err = norm(res);
    if show
        fprintf('step %3d | err: %8.3e | lambda (%16.9e,%16.9e) | mu (%16.9e,%16.9e) \n',step,err,real(k),imag(k),real(mu),imag(mu));
    end
    if err<reltol
        isConverged = true;
        break
    end
    upd = -J\res;
    x = x + upd(1:n);
    y = y + upd(n+1:2*n);
    k = k + upd(2*n+1);
    mu = mu + upd(2*n+2);
    x = x/norm(x);
    y = y/norm(y);
    y = conj(y);
end

% return values:
w = sqrt(mu); % angular frequency
if ~isConverged || k < kmin % return nan if not converge or cutoff was found
    k = nan;
    w = nan;
    x = nan(size(x));
    y = nan(size(y));
end
