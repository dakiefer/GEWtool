classdef Geometry
% Geometry - Describes a 1D layered geometry (mesh).
% Collects all quantities describing the geometry and mesh. Maps local degrees of 
% freedom to global ones. Remembers boundary condition setup. There is usually no 
% need to use this class explicitly (used internally by Waveguide).
% This could probably be replaced by a standard finite element mesh.
%
% See also Plate, Cylinder, Waveguide.
% 
% 2022 - Daniel A. Kiefer
% Institut Langevin, Paris, France
% 
    
properties 
    nLay         % number of layers
    nItf         % number of interfaces, includes outer ones
    N            % number of nodal points for each layer
    Nudof        % number of displacement dofs (e.g., ux, uy, uz) for each layer
    yItf         % position of interfaces in meter [nLay x 2]
    y            % cell array with nodal points for each layer in meter
    h            % thickness for each layer in meter
    % nodesOfElem  % connectivity map: row e contains left and right node num of elem e -> TODO not being used
    ldofBC       % cell array of local dofs at boundaries: e.g. [1, N; N+1, 2*N] -> [upper, lower]
    gdofDBC      % fixed degrees of freedom
    gdofOfLay    % cell array with global dofs assigned to layers
    Ndof         % total number of degrees of freedom
end

properties (Dependent)
    gdofBC      % a link to all gdofOfLay(ldofBC) assembled into one array: size [ui, 2 (bottom/top)]
end

methods

    function obj = Geometry(yItf, N, Nudof)
        % Geometry - Create a Geometry object for a 1D layered structure (mesh).
        % Arguments:
        % - yItf:  coordinates of interface in meter
        %          Can be a vector or of size [Nlay x 2] (lower/upper limit for 
        %          every layer).
        % - N:     discretization order for each layer [1 x Nlay]
        % - Nudof: displacement digrees of freedom for each layer [1 x Nlay].
        %          (optional, default per layer: 3)

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
        for e=1:length(N) % length(N) == number of layers
            obj.ldofBC{e} = (0:Nudof(e)-1).'*N(e) + [1, N(e)];
            etad = Layer.nodes(N(e)); 
            obj.y{e} = etad*obj.h(e) + obj.yItf(e, 1);
            if e == 1
                gdofE = 1:Nudof(e)*N(e);
            else
                Ncommon = min(Nudof(e-1:e)); % number of common degrees of freedom at interface
                gdofBCLast = obj.gdofOfLay{e-1}(obj.ldofBC{e-1}); % global dofs at boundary of last layer
                gdofCommon = gdofBCLast(end-Ncommon+1:end,2).'; % dofs of previous layer at interface
                gdofLastLay = obj.gdofOfLay{e-1}(end);
                Nnew = Nudof(e)*N(e)-Ncommon; % subtract common number of dofs with previous layer
                gdofNew = gdofLastLay + (1:Nnew); 
                ldofInd = 1:Nudof(e)*N(e);
                ldofCommon = obj.ldofBC{e}(1:Ncommon);
                ldofNew = setdiff(ldofInd, ldofCommon);
                gdofE = zeros(1,length(ldofInd));
                gdofE(ldofNew) = gdofNew; gdofE(ldofCommon) = gdofCommon;
            end
            obj.gdofOfLay{e} = gdofE;
        end
        obj.Ndof = numel(unique([obj.gdofOfLay{:}]));
    end

    function gdofBC = get.gdofBC(obj)
        gdofBC = cell(length(obj.N), 1); % initialize
        for e=1:length(obj.N) % loop over layers 
            gdofE = obj.gdofOfLay{e};
            gdofBC{e} = gdofE(obj.ldofBC{e});
        end
    end

end % methods

end % class
