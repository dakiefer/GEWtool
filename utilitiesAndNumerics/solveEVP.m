function [lbd, y] = solveEVP(H, nModes, opts)
    if nargout == 1
        if opts.subspace
            lbd = eigs(H, nModes, opts.target);
        else
            lbd = eig(H, 'vector');
        end
    else
        if opts.subspace
            [y, lbd] = eigs(H, nModes, opts.target);
            lbd = diag(lbd);
        else
            [y, lbd] = eig(H, 'vector');
        end
    end
end
