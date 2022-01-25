classdef Geometry
    
properties 
    nElem
    nNodes
    N
    Nudof
    yItf
    y
    h
    nodesOfElem
    ldofBC
    gdofOfElem
end

methods

    function geom = Geometry(yItf, N, Nudof)
        % mesh parameters:
        if nargin < 3, Nudof = 3*ones(size(N)); end
        geom.nElem = length(yItf) - 1;
        geom.N = N(:); % discretization order: number of collocation points
        geom.Nudof = Nudof(:); % how many of the ux, uy, uz degrees of freedom
        geom.nNodes = length(yItf);
        geom.yItf = [yItf(1:end-1).', yItf(2:end).'];
        geom.h = diff(geom.yItf, 1, 2); % first order along 2nd dim
        % mesh connectivity:
        geom.nodesOfElem = zeros(geom.nElem, 2); % connectivity map
        geom.nodesOfElem(:,1) = 1:geom.nNodes-1;
        geom.nodesOfElem(:,2) = 2:geom.nNodes;
        % boundary condition nodes:
        dof0 = 0;
        for e=1:length(N)
            geom.ldofBC{e} = (0:Nudof(e)-1)*N(e) + [1; N(e)];
            [etad, ~] = chebdif(N(e), 1);
            geom.y{e} = (geom.h(e)*etad + geom.yItf(e, 2) + geom.yItf(e, 1))/2;
            geom.gdofOfElem{e} = dof0 + (1:Nudof(e)*N(e)); dof0 = dof0 + Nudof(e)*N(e);
        end
    end

end % methods

end % class
