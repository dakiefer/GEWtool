classdef PlateLeaky < Plate
% PlateLeaky - Represents quasi-guided waves in a plate loaded by a half-space.
% Displacement ansatz: u(x,y,z,t) = u(z)*exp(i k x - i w t)
% 
% Example:
% mat = Material('steel'); % load material data
% fluid.c = 1500; fluid.rho = 1000;
% h = 1e-3; % thickness in m
% N = 20; % discretization (number of nodal points)
% plate = PlateLeaky({mat fluid}, [h inf], N); % create waveguide description
% 
% See also Plate.Plate, Cylinder, Waveguide.
% 
% 2025 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France

properties
    halfSpaces = [] % TODO field "loading" as a array of struct with load.mat and load.at = "top"/"bottom"
    % TODO expand "loading" struct array in terms of wave velocities (2 for solid loading) 
end

methods
	function obj = PlateLeaky(mats, zs, Ns)
        if ~iscell(mats)
            mats = num2cell(mats);
        end
        halfSpaces = PlateLeaky.parseLoading(mats,zs);
        mats = mats{~isinf(zs)}; zs = zs(~isinf(zs)); % crop to finite layers
        obj = obj@Plate(mats, zs, Ns);
        obj.halfSpaces = halfSpaces;
    end
    function obj = assembleLayers(obj, udof, n)
        obj = assembleLayers@Waveguide(obj, udof, n);
        for i = 1:length(obj.halfSpaces)
            loading = obj.halfSpaces(i);
            Ndof = size(obj.op.L0,1); % increases with every iteration
            [dofA, dofU] = PlateLeaky.getCouplingDOFs(loading,obj.geom,Ndof);
            obj = incorporateLoading(obj, loading, dofA, dofU); 
        end
        if length(obj.halfSpaces) == 2 && obj.halfSpaces(1).mat == obj.halfSpaces(2).mat
            obj.op.Rtop = obj.op.Rtop - obj.op.Rbottom; % top - bottom (waves radiated away from the plate)
            obj.op = rmfield(obj.op,'Rbottom'); 
        end
    end
    function obj = incorporateLoading(obj, loading, dofA, dofU)
        if loading.at == "top" % different signs at top and bottom
            sig = 1; 
        else
            sig = -1; 
        end
        % sides = {'a', 'b'}; % used for labeling the sides
        op = obj.op;
        nDof = max(dofA); % new total number of DOFs

        % expand matrices: 
        opNameList = fieldnames(op);
        for i=1:length(opNameList)
            op.(opNameList{i})(nDof,nDof) = 0; 
        end
        R = zeros(nDof);  % radiation matrix (models nonpolynomial terms)
        
        if isa(loading.mat,'MaterialFluid')
            % continuity of normal displacements: ibeta*A - uz = 0
            op.L0(dofA,dofU) = -1;
            R(dofA,dofA)     = +1;
            % balance of tractions: add the boundary term [v*tA] to the FE matrices. 
            % the traction induced by the fluid is tA = -w^2*rhoA*ez*A
            op.M(dofU,dofA) = sig*loading.mat.rho/obj.np.rho0; % normalized mass density
            op.("R"+loading.at) = R; 
        elseif isa(loading.mat,'MaterialIsotropic')
            error('GEWTOOL:not implemented yet.') % TODO implement
        end
        obj.op = op; 
    end
    function op = opExpandTerm(obj, opName, extMat)
        op = obj.op;
        L2 = op.L2; L1 = op.L1; L0 = op.L0; M = op.M; R = op.(opName);
        op = rmfield(op, {'L2','L1','L0','M',char(opName)}); 
        Z = zeros(size(L0));
        cf = extMat.cl/obj.np.fh0; % normalized wave speed
        LL2 = [L2  ,   -R   ;
               Z   ,   L2    ];
        LL1 = blkdiag(L1,L1);
        LL0 = [L0  ,   Z    ; 
               R   ,   L0   ];
        MM  = [M   ,   -1/cf^2*R  ; 
               Z   ,      M       ];
        % expand remaining R-matrices: 
        opList = fieldnames(op);
        for i=1:length(opList)
            Ri = op.(opList{i});
            op.(opList{i}) = blkdiag(Ri, Ri); % expand
        end
        op.L2 = LL2; op.L1 = LL1; op.L0 = LL0; op.M = MM;
    end
    function gew = symmetrizeGeometry(obj)
        % symmetrizeGeometry - upper symmetric half of the original plate.
        % Creates a new Plate object describing only the upper symmetric half 
        % of the original Plate object. 
        % Throws a warning if the plate is not symmetric in geometry and materials. 
        % This function is used by LambS, LambA, LambSA, etc., and you will 
        % usually not need to call it explicitly.
        if ~obj.decouplesSA('v')
            error('GEWTOOL:symmetrizeGeometry','The setup is not symmetric. I cannot symmetrize the geometry.'); 
        end
        gew = symmetrizeGeometry@Plate(obj);
        gew.exteriorMat = obj.exteriorMat;
        gew.exteriorMat{1} = []; % loading only at top side
    end
    function decoupl = decouplesSA(obj, verb)
        % decouplesSA - Tests whether symmetric and antisymmetric waves decouple.
        % Usage: 
        % decoupl = decouplesSA;       Returns true if SA waves decouple.
        % decoupl = decouplesSA('v');  Throw warning indicating reason.
        % 
        % See also: decouplesLambvsSH.

        % verify symmetry of materials:
        if nargin ~= 2
            verb = 'nonVerb';
        end
        bothSides = ~isempty(obj.exteriorMat{1}) && ~isempty(obj.exteriorMat{2});
        bothSame = bothSides && obj.exteriorMat{1} == obj.exteriorMat{2}; % short circuit to avoid error for empty entries
        decoupl = bothSame && decouplesSA@Plate(obj,verb);
    end
end

methods (Static)
    function [dofA, dofU] = getCouplingDOFs(halfspace,geom,Ndof)
        if halfspace.at == "top"
            lay = geom.nLay; side = 2;
        elseif halfspace.at == "bottom" 
            lay = 1; side = 1;
        end
        if isa(halfspace.mat,'MaterialFluid') % is a fluid
            dofU = geom.gdofBC{lay}(end,side); % end -> last displacement component is always normal to the plate
            dofA = Ndof + 1; 
        elseif isa(halfspace.mat,'MaterialIsotropic')
            dofU = geom.gdofBC{lay}(:,side);
            dofA = Ndof + (1:length(dofU)); % same number of additional unknows as number of displacement components
        end
    end
    function halfSpaces = parseLoading(matList,zs)
        ind = sort(find(isinf(zs))); 
        if length(ind) > 2
            error('GEWTOOL:PlateLeaky','The plate cannot be loaded with more than two halfspaces.');
        elseif isempty(ind)
            error('GEWTOOL:PlateLeaky','Provide at least one loading halfspace or use the Plate class for nonleaky waves.');
        end
        for i = 1:length(ind)
            mati = matList{ind(i)};
            if ~isa(mati,'MaterialIsotropic') && ~isa(mati,'MaterialFluid') % only these are supported for now
                warning('GEWTOOL:parseLoading','Trying to convert loading Material to class "MaterialIsotropic". To hide this warning, load your material with "MaterialIsotropic". Loading with anisotropic materials is not supported.');
                mati = MaterialIsotropic(mati); % convert from "Material" to "MaterialIsotropic"
            end
            halfSpaces(i).mat = mati;
            if ind(i) == 1
                halfSpaces(i).at  = "bottom";
            else
                halfSpaces(i).at = "top";
            end
        end
    end
end

end