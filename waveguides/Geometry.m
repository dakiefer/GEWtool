classdef Geometry
    
properties 
    nLay         % number of layers
    nItf         % number of interfaces, includes outer ones
    N            % number of collocation points for each layer
    Nudof        % number of displacement dofs for each layer
    yItf         % position of interfaces in m
    y            % cell array with collocation points for each layer in m 
    h            % thickness for each layer in m
    % nodesOfElem  % connectivity map: row e contains left and right node num of elem e -> TODO not being used
    ldofBC       % local dofs of BC: e.g. [1, N+1; N, 2*N] -> [upper; lower]
    gdofOfLay    % cell array with global dofs assigned to layers
end

methods

    function geom = Geometry(yItf, N, Nudof)
        % mesh parameters:
        if nargin < 3, Nudof = 3*ones(size(N)); end
        geom.nLay = length(yItf) - 1;
        geom.N = N(:); % discretization order: number of collocation points
        geom.Nudof = Nudof(:); % how many of the ux, uy, uz degrees of freedom
        geom.nItf = length(yItf);
        geom.yItf = [yItf(1:end-1).', yItf(2:end).'];
        geom.h = diff(geom.yItf, 1, 2); % first order along 2nd dim
        % mesh connectivity:
        % geom.nodesOfElem = zeros(geom.nLay, 2); % connectivity map
        % geom.nodesOfElem(:,1) = 1:geom.nItf-1;
        % geom.nodesOfElem(:,2) = 2:geom.nItf;
        % boundary condition nodes:
        dof0 = 0;
        for e=1:length(N)
            geom.ldofBC{e} = (0:Nudof(e)-1)*N(e) + [1; N(e)];
            [etad, ~] = chebdif(N(e), 1);
            geom.y{e} = (geom.h(e)*etad + geom.yItf(e, 2) + geom.yItf(e, 1))/2;
            geom.gdofOfLay{e} = dof0 + (1:Nudof(e)*N(e)); dof0 = dof0 + Nudof(e)*N(e);
        end
    end

end % methods

end % class
