function opts = parseSolverOpts(opts, nModes)
% DEFAULTSOLVEROPTS - Parse the solver options.
% Returns a structure specifying the options for the solvers computeK() and
% computeW(). The user specified options are extended by the default values of
% the non-specified options.
% 
% 2022-2023 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France

if ~isfield(opts, 'sparse'),   opts.sparse = false;    end
if ~isfield(opts, 'subspace'), opts.subspace = false;  end
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
        opts.target = "smallestabs";
    end
end
if ~isempty(nModes) && ~isinf(nModes) && ( ~isscalar(nModes) || nModes ~= round(nModes) )
    error('GEWTOOL:wrongArg', 'Argument "nModes" must be a scalar integer.');
end

end