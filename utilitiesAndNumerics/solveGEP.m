function [lbd, eVec] = solveGEP(A, B, nModes, opts)
    if nargout == 1
        if opts.subspace
            lbd = eigs(A, B, nModes, opts.target);
        else
            lbd = eig(A, B, 'chol', 'vector');
        end
    else
        if opts.subspace
            [eVec, lbd] = eigs(A, B, nModes, opts.target);
            lbd = diag(lbd);
        else
            [eVec, lbd] = eig(A, B, 'chol', 'vector');
        end
    end
end