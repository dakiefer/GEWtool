function [opts, nModes] = parseSolverOpts(opts, op, nModes)
% DEFAULTSOLVEROPTS - Parse the solver options.
% Returns a structure specifying the options for the solvers computeK() and
% computeW(). The user specified options are extended by the default values of
% the non-specified options.
% 
% 2022-2023 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France

if ~isfield(opts, 'show'),     opts.show = false;      end
if opts.compute == "k" 
    isSingular = rcond(op.L0 + op.M) < 1e2*eps; % k from -[L0 + w^2*M]*u = [ (ik)^2*L2 + ik*L1 ]*u, test at w = 1
elseif opts.compute == "w"
    isSingular = rcond(op.L2 + op.L1 + op.L0) < 1e2*eps; % w from: -[(ik)^2*L2 + ik*L1 + L0]*u = w^2*M*u, test at k = 1
end
if ~isfield(opts, 'subspace') % guess best algorithm to use
    if ~isSingular && size(op.M,1) > 60 % 60 was determined empirically
        opts.subspace = true;  % default for big matrices
    else
        opts.subspace = false; % default for small matrices
    end
else
    if opts.subspace && isSingular && ( ~isfield(opts, 'target') || ~isnumeric(opts.target) )
        warning('GEWTOOL:singular', 'The eigenvalue problem is singular or close to singular. If you encounter stability issues, switch to the QZ-algorithm by providing opts.subspace = false to computeK(), computeW(). Alternatively, try providing a target in opts.target. I will proceed anyways.');
    end
end
if ~isfield(opts, 'sparse'), opts.sparse = opts.subspace; end % default: sparse when subspace methods are used
if (isempty(nModes) || isinf(nModes)) && opts.subspace % default for subspace
    if size(op.M,1) <= 20
        nModes = size(op.M,1);
    else 
        nModes = 20; 
        warning('GEWTOOL:missingNModes', 'You should provide the number of modes you need, e.g., computeW(gew, k, nModes). I will be computing %d modes per default.', nModes);
    end
elseif (isempty(nModes) || isinf(nModes)) && ~opts.subspace % default for QZ
    nModes = size(op.M,1); 
end
if nModes > size(op.M,1) % if the user explicitly provided nModes bigger than the matrix
    warning('GEWTOOL:computeW:tooManyModes', 'More modes requested than available. Resetting nModes to the matrix size.')
    nModes = size(op.M,1);
end
if ~isfield(opts, 'parallel')  
    if ~canParallel() % check first if parallel computing toolbox is available
        opts.parallel = false; 
    else
        opts.parallel = ~isempty(gcp('nocreate')); % if a parallel pool exists, then use it!
    end
end
if islogical(opts.parallel) % convert logical value to number of workers
    if opts.parallel, opts.parallel = inf; else, opts.parallel = 0; 
end 
if ~opts.subspace && opts.sparse
    warning('GEWTOOL:ignoringSparse',...
        'Sparse matrices are only supported in combination with the subspace solver, i.e., eigs(). Switching to subspace method. To hide this message set opts.subspace=true;');
    opts.subspace = true; % switch to eigs()
end 
if opts.subspace 
    if ~isfield(opts,'target')
        opts.target = 'smallestabs';
    end
end
if ~isempty(nModes) && ~isinf(nModes) && ( ~isscalar(nModes) || nModes ~= round(nModes) )
    error('GEWTOOL:wrongArg', 'Argument "nModes" must be a scalar integer.');
end
if ~isfield(opts, 'eigenvecs'), opts.eigenvecs = true; end % default: compute eigenvectors
if ~isfield(opts, 'standardEVP'), opts.standardEVP = isdiag(op.M) && rcond(op.M) >= 1e-4; end % default value
if ~isfield(opts, 'trace'), opts.trace = true; end % default: reorder modes so that they are properly sorted (grouped)

if opts.show, disp(opts); end % display options for debugging

end