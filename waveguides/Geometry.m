classdef Geometry
    
properties 
    nLay         % number of layers
    nItf         % number of interfaces, includes outer ones
    N            % number of collocation points for each layer
    Nudof        % number of displacement dofs for each layer
    yItf         % position of interfaces in m [nLay x 2]
    y            % cell array with collocation points for each layer in m 
    h            % thickness for each layer in m
    % nodesOfElem  % connectivity map: row e contains left and right node num of elem e -> TODO not being used
    ldofBC       % local dofs of BC: e.g. [1, N+1; N, 2*N] -> [upper; lower]
    gdofBC = []; % global dofs of BCs: differs from ldofBC for multilayer problems
    gdofOfLay    % cell array with global dofs assigned to layers
end

methods

    function obj = Geometry(yItf, N, Nudof)
        % mesh parameters:
        if nargin < 3, Nudof = 3*ones(size(N)); end
        if ~isvector(yItf)
            if size(yItf, 2)~=2
                error('Incorrect size of interface coordinates. Expected a vector or nx2 array.');
            end
            yItf = [yItf(:,1).', yItf(end,2)];
        end
        obj.nLay = length(yItf) - 1;
        obj.N = N(:); % discretization order: number of collocation points
        obj.Nudof = Nudof(:); % how many of the ux, uy, uz degrees of freedom
        obj.nItf = length(yItf);
        obj.yItf = [yItf(1:end-1).', yItf(2:end).'];
        obj.h = diff(obj.yItf, 1, 2); % first order along 2nd dim
        % mesh connectivity:
        % obj.nodesOfElem = zeros(obj.nLay, 2); % connectivity map
        % obj.nodesOfElem(:,1) = 1:obj.nItf-1;
        % obj.nodesOfElem(:,2) = 2:obj.nItf;
        % boundary condition nodes:
        dof0 = 0;
        for e=1:length(N) % length(N) == number of layers
            obj.ldofBC{e} = (0:Nudof(e)-1)*N(e) + [1; N(e)];
            [etad, ~] = chebdif(N(e), 1);
            obj.y{e} = (-obj.h(e)*etad + obj.yItf(e, 2) + obj.yItf(e, 1))/2;
            obj.gdofOfLay{e} = dof0 + (1:Nudof(e)*N(e)); dof0 = dof0 + Nudof(e)*N(e);
            obj.gdofBC = [obj.gdofBC, obj.gdofOfLay{e}(obj.ldofBC{e})];
        end
    end

end % methods

end % class
