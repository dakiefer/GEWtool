function [u, k] = GEWpolyeig(P0, P1, P2, nModes, opts)
% GEWpolyeig - Solve the quadratic eigenvalue problem [k^2 P2 + k P1 + P0]*u = 0.
% GEWpolyeig is a bit faster than polyeig() and supports subspace methods as
% well as sparse matrices.
%
% Arguments:
% - P0, P1, P2:   Coefficient matrices of the quadratic eigenvalue problem
% - nModes:  (optional) Number of eigensolutions (smallest eigenvalue magnitude).
%
% Return values:
% - u:  The eigenvectors as a matrix of size [size(P0,1), nModes]
% - k:  The eigenvalues as a vector of size [nModes, 1]
% 
% See also polyeig.
% 
% 2023 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France

if nargin < 5, opts = []; end
if ~isfield(opts, 'sparse'),   opts.sparse = false;    end
if ~isfield(opts, 'subspace'), opts.subspace = false;  end
if ~isfield(opts, 'parallel'), opts.parallel = false;  end
if ~opts.subspace && opts.sparse
    warning('GEWTOOL:ignoringSparse',...
        'Sparse matrices are only supported in combination with the subspace solver, i.e., eigs(). Switching to subspace method. To hide this message set opts.subspace=true;');
    opts.subspace = true; % switch to eigs()
end 
if opts.parallel, opts.parallel = inf; else, opts.parallel = 0; end % set the number of workers
if nargin < 4, nModes = size(P0,1); end

% N = size(varargin{1}); % size of matrices
% p = nargin-1; % polynomial order

% % linearization preserving structure:
% % % quadratic only for now, i.e., polynomial: lbd^2*A + lbd*B + C
% % linearize as A*x = lbd*B*x
% r = 1; % any real value should work
% B = [1i*P2,   -r*P2; 
%      r*P2,    r*P1+1i*P0 ];
% A= -[r*P2 + 1i*P1,  1i*P0;
%       -1i*P0,       r*P0  ];

% % companion matrix linearization
p = 2; % polynomial order (only two for now)
n = size(P0,1); 

% allocate matrices:
if opts.sparse
    A = speye(n*p);
    B = sparse(n*p, n*p); % sparse zeros-matrix
else
    A = eye(n*p);
    B = zeros(n*p);
end

A(1:n,1:n) = P0;         % P0 in first block
B(n+1:n*p+1:end) = 1;    % I on lower block-diagonal
B(1:n,1:n) = -P1;        % -P1 in block (1,1)
B(1:n,n+1:2*n) = -P2;    % -P2 in block (1,2)

if opts.subspace
    [up, k] = eigs(A, B, nModes, "smallestabs"); % complex wavenumbers -> non-Hermitian
    k = diag(k);
else
    [up, k] = eig(A, B, 'vector'); % complex wavenumbers -> non-Hermitian
end
u = up(1:n,:); % crop to first part of eigenvector.
% NOTE: The cropping is done differently in polyeig()! The cropping done here
% passes the test in tests/test_GEWdat_plate_field_calculation.m

end
