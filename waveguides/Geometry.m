classdef Geometry
% Geometry - Describes a 1D layered geometry (mesh).
% Collects all quantities describing the geometry and mesh. Maps local degrees of 
% freedom to global ones. Remembers boundary condition setup. There is usually no 
% need to use this class explicitly (used internally by Waveguide).
% This could probably be replaced by a standard finite element mesh.
%
% See also Plate, Cylinder, Waveguide.
% 
% 2022 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France
    
properties 
    N            % number of nodal points for each layer
    Nudof        % number of degrees of freedom (displacement components, elec. potential, etc...)
    zItf         % position of interfaces in meter [nLay x 2]
    z            % cell array with nodal points for each layer in meter
    % nodesOfElem  % connectivity map: row e contains left and right node num of elem e -> TODO not being used
    gdofDBC      % fixed degrees of freedom
    symmetrized = false  % true/false: whether the geometry represents the full or symmetry part
end

properties (Dependent)
    nLay         % number of layers
    nItf         % number of interfaces, includes outer ones
    hl           % thickness for each layer in meter
    Ndof         % total number of degrees of freedom
    gdofOfLay    % cell array with global dofs assigned to layers
    ldofBC       % cell array of local dofs at boundaries: e.g. [1, N; N+1, 2*N] -> [upper, lower]
    gdofBC       % a link to all gdofOfLay(ldofBC) assembled into one array: size [ui, 2 (bottom/top)]
    gdofFree     % degrees of freedom where no Dirichlet-BCs are imposed
    gdofRedX     % degrees of freedom for the x-displacements in the reduced matrices (due to Dirichlet BC)
    gdofRedZ     % degrees of freedom for the z-displacements in the reduced matrices (due to Dirichlet BC)
end

methods

    function obj = Geometry(zItf, N, Nudof)
        % Geometry - Create a Geometry object for a 1D layered structure (mesh).
        % Arguments:
        % - zItf:  coordinates of interface in meter
        %          Can be a vector or of size [Nlay x 2] (lower/upper limit for 
        %          every layer).
        % - N:     discretization order for each layer [1 x Nlay]
        % - Nudof: displacement digrees of freedom for each layer [1 x Nlay].
        %          (optional, default per layer: 3)

        % mesh parameters:
        if nargin < 3, Nudof = []; end
        if ~isvector(zItf) % nLay x 2 given
            if size(zItf, 2)~=2
                error('Incorrect size of interface coordinates. Expected a vector or nx2 array.');
            end
            zItf = [zItf(:,1).', zItf(end,2)]; % reduce to vector nLay+1 x 1
        end
        obj.zItf = [zItf(1:end-1).', zItf(2:end).'];
        obj.N = N(:); % discretization order: number of collocation points
        obj.Nudof = Nudof(:); % how many of the ux, uy, uz degrees of freedom
        for e=1:obj.nLay % length(N) == number of layers
            etad = Layer.nodes(obj.N(e)); 
            obj.z{e} = etad*obj.hl(e) + obj.zItf(e, 1);
        end
        % mesh connectivity:
        % obj.nodesOfElem = zeros(obj.nLay, 2); % connectivity map
        % obj.nodesOfElem(:,1) = 1:obj.nItf-1;
        % obj.nodesOfElem(:,2) = 2:obj.nItf;
    end

    function nLay = get.nLay(obj)
        nLay = size(obj.zItf,1);
    end

    function nItf  = get.nItf(obj)
        nItf = size(obj.zItf,1)+1;
    end

    function hl = get.hl(obj)
        hl = diff(obj.zItf, 1, 2); % first order along 2nd dim
    end

    function Ndof = get.Ndof(obj)
        Ndof = numel(unique([obj.gdofOfLay{:}]));
    end

    function gdofOfLay = get.gdofOfLay(obj)
        % boundary condition nodes:
        gdofOfLay = cell(obj.nLay,1);
        for e=1:obj.nLay
            if e == 1
                gdofE = 1:obj.Nudof(e)*obj.N(e);
            else
                Ncommon = min(obj.Nudof(e-1:e)); % number of common degrees of freedom at interface
                gdofBCLast = gdofOfLay{e-1}(obj.ldofBC{e-1}); % global dofs at boundary of last layer
                gdofCommon = gdofBCLast(end-Ncommon+1:end,2).'; % dofs of previous layer at interface
                gdofLastLay = gdofOfLay{e-1}(end);
                Nnew = obj.Nudof(e)*obj.N(e)-Ncommon; % subtract common number of dofs with previous layer
                gdofNew = gdofLastLay + (1:Nnew); 
                ldofInd = 1:obj.Nudof(e)*obj.N(e);
                ldofCommon = obj.ldofBC{e}(1:Ncommon);
                ldofNew = setdiff(ldofInd, ldofCommon);
                gdofE = zeros(1,length(ldofInd));
                gdofE(ldofNew) = gdofNew; gdofE(ldofCommon) = gdofCommon;
            end
            gdofOfLay{e} = gdofE;
        end
    end

    function ldofBC = get.ldofBC(obj)
        ldofBC = cell(obj.nLay,1);
        for e=1:obj.nLay
            ldofBC{e} = (0:obj.Nudof(e)-1).'*obj.N(e) + [1, obj.N(e)];
        end
    end

    function gdofBC = get.gdofBC(obj)
        gdofBC = cell(length(obj.N), 1); % initialize
        for e=1:length(obj.N) % loop over layers 
            gdofE = obj.gdofOfLay{e};
            gdofBC{e} = gdofE(obj.ldofBC{e});
        end
    end

    function gdofFree = get.gdofFree(obj)
        dofs = [obj.gdofOfLay{:}];
        gdofFree = setdiff(dofs, obj.gdofDBC);
    end

    function dofx = get.gdofRedX(obj)
        if obj.nLay ~= 1; error('GEWTOOL:Geometry','Not implemented for Multilayers.'); end
        if obj.Nudof == 1 % SH waves are modeled
            error('GEWTOOL:Geometry','There are no ux-displacements to extract.');
        end
        dofAll = 1:length(obj.gdofFree);
        indx = obj.gdofFree <= obj.N;
        dofx = dofAll(indx);
    end

    function dofz = get.gdofRedZ(obj)
        if obj.nLay ~= 1; error('GEWTOOL:Geometry','Not implemented for Multilayers.'); end
        if obj.Nudof == 2 % Lamb waves 
            dofMin = obj.N+1; dofMax = 2*obj.N; 
        elseif obj.Nudof >= 3 % fully-coupled waves or piezoelectric waves
            dofMin = 2*obj.N+1; dofMax = 3*obj.N; 
        else
            error('GEWTOOL:Geometry','There are no uz-displacements to extract.');
        end
        dofAll = 1:length(obj.gdofFree);
        indz = obj.gdofFree >= dofMin & obj.gdofFree <= dofMax;
        dofz = dofAll(indz);
    end

    function [l, n, u] = layerNodeCompOf(obj, gdof)
    % layerNodeCompOf - find the layer number l and node number n corresponding to
    % the global degree of freedom gdof. 
        gdof = obj.gdofFree(gdof); % get gdof before removal of Dirichlet nodes 
        for l = 1:obj.nLay
            gdofList = obj.gdofOfLay{l}; 
            nxu = find(gdofList == gdof);
            if ~isempty(nxu)
                n = mod(nxu-1,obj.N)+1; 
                u = floor((nxu-1)/obj.N)+1;
                return;
            end
        end
    end

end % methods

end % class
