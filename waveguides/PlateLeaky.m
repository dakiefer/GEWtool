classdef PlateLeaky < PlateClosed
% PlateLeaky - Represents quasi-guided waves in a plate loaded by a half-space.
% 
% Quasi-guided waves encompass "leaky waves" and "trapped waves" (= surface
% waves or quasi-Scholte waves). 
% 
% Usage:
% -----------------------------------------------
% mat = Material('steel');          % load material data
% fluid.c = 1500; fluid.rho = 1000; % wave speed and density
% h = 1e-3; % thickness in m
% N = 20;   % discretization (number of nodal points)
% plate = Plate({mat fluid}, [h inf], N); % waveguide description
% -----------------------------------------------
% 
% The displacement ansatz for quasi-guided waves in the plate is
% 
% u(x,y,z,t) = u(z)*exp(i k x - i w t)
% 
% and in the exterior half-spaces: 
% 
% u(x,y,z,t) = ∑j Aj exp(i betaj z)*exp(i k x - i w t).
% 
% Quasi-guided waves are governed by the following EVP that is nonlinear in the
% eigenvalue ik: 
% 
% [ (ik)^2*L2 + ik*L1 + L0 + w^2*M + ∑j ibetaj*Rj ]*q = 0,
% 
% where q:       eigenvector containing the displacements in the plate and the "nA"
%                bulk wave amplitudes "Aj".
%       ik:      eigenvalue -> horizontal wavenumber (times i)  
%       ibetaj:  vertical wavenumber (times i) of the jth bulk wave
%       w:       angular frequency (parameter) 
%       Li,M,Rj: n x n-matrices stored in "obj.gew.opNonlin".
% 
% The above problem is nonlinear in ik because the ibetaj are related to ik by
% 
% ik^2 + ibetaj^2 = (iw/cj)^2 with bulk wave velocities "cj". 
% 
% The nonlinear eigenvalue problem is transformed to a polynomial one by
% introducing the new eigenvectors:
% 
% Psi = [ ibetaj q ]
%       [    q     ] 
% 
% The above is applied recursively for each bulk wave j, i.e., with 2 bulk waves
% the eigenvector has four blocks: 
% 
% Psi = [ ibeta1 ibeta2 q ]  <- block 1
%       [    ibeta2 q     ]  <- block 2
%       [    ibeta1 q     ]  <- block 3
%       [        q        ]  <- block 4
% 
% For radiation into fluid media, the above transformation leads to a quadratic
% eigenvalue problem while radiation into solids leads to a cubic eigenvalue
% problem. In the latter case, GEWtool applies a partial companion linearization
% to reduce the problem to a quadratic one. In both cases, GEWtool finaly solves
% an eigenvalue problem of the form
% 
% [ (ik)^2*L2 + ik*L1 + L0 + w^2*M ]*Psi = 0
% 
% with Psi: eigenvector stored in "obj.Psi"
%      ik: eigenvalue (as before) stored in "obj.k"
%      Li,M: matrices stored in "obj.gew.op"
%
% See also PlateLeaky.PlateLeaky, Cylinder, Waveguide.
% 
% 2025 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France

properties
    halfSpaces = [] % TODO field "loading" as a array of struct with load.mat and load.at = "top"/"bottom"
    opNonlin = [] % matrices of the nonlinear eigenvalue problem
    % TODO expand "loading" struct array in terms of wave velocities (2 for solid loading) 
end

methods
	function obj = PlateLeaky(mats, zs, Ns, halfSpaces)
        obj = obj@PlateClosed(mats, zs, Ns);
        for i = 1:length(halfSpaces) 
            halfSpaces(i).eliminated = false; % augument this field to remember if the half space has already been eliminated
        end
        obj.halfSpaces = halfSpaces;
    end
    function obj = assembleLayers(obj, udof, n)
        obj = assembleLayers@Waveguide(obj, udof, n);
        obj.opNonlin = obj.op;
        % setup radiation matrices of the nonlinear eigenvalue problem:
        for i = 1:length(obj.halfSpaces)
            loading = obj.halfSpaces(i);
            Ndof = size(obj.opNonlin.L0,1); % increases with every iteration
            [dofA, dofU] = PlateLeaky.getCouplingDOFs(loading, obj.geom, Ndof);
            obj.halfSpaces(i).dofA = dofA; % save for postprocessing purposes 
            obj.halfSpaces(i).dofU = dofU; % save for postprocessing purposes 
            obj = incorporateLoading(obj, loading, dofA, dofU, udof); 
        end 
        % simplify to leaky-only when the halfspaces on both sides are equal: 
        if length(obj.halfSpaces) == 2 && obj.halfSpaces(1).mat == obj.halfSpaces(2).mat
            obj.opNonlin.Rtop = obj.opNonlin.Rtop - obj.opNonlin.Rbottom; % top - bottom (waves radiated away from the plate)
            obj.opNonlin = rmfield(obj.opNonlin,'Rbottom'); 
            obj.halfSpaces(1).eliminated = true; %  = obj.halfSpaces(2); 
        end
        % transform into a polynomial eigenvalue problem in a higher-dimensional
        % state space: 
        op = obj.opNonlin; 
        for i = 1:length(obj.halfSpaces)
            op = getPolynomialForm(obj, op, obj.halfSpaces(i));
        end
        op = PlateLeaky.reduceToQuadratic(op); 
        obj.op = op;
    end
    function obj = incorporateLoading(obj, loading, dofA, dofU, udof)
        % sides = {'a', 'b'}; % used for labeling the sides
        nDof = max(dofA); % new total number of DOFs

        % expand matrices: 
        opNameList = fieldnames(obj.opNonlin);
        for i=1:length(opNameList)
            obj.opNonlin.(opNameList{i})(nDof,nDof) = 0; 
        end
        
        % build radiation matrices:
        if isa(loading.mat,'MaterialFluid')
            obj = incorporateFluidLoading(obj, loading, dofA, dofU);
        elseif isa(loading.mat,'MaterialIsotropic')
            obj = incorporateSolidLoading(obj, loading, dofA, dofU, udof); 
        end
    end
    function obj = incorporateFluidLoading(obj, loading, dofA, dofU)
        % given degrees of freedom (dofA: additional for halfspace, dofU: displacements at boundary)
        opN = obj.opNonlin;
        nDof = size(opN.M,1);
        % allocate new matrix:
        R = zeros(nDof);  % radiation matrix (models nonpolynomial terms)
        % get coupling matrices: 
        coupl = PlateLeaky.couplingMatricesFluid(loading,obj.np);

        % continuity of normal displacements: ibeta*A - uz = 0
        opN.L0(dofA,dofU) = -1;
        R(dofA,dofA)      = coupl.Ug; % = 1
        
        % balance of tractions: add the boundary term [v*tA] to the FE matrices. 
        % the traction induced by the fluid is tA = -w^2*rhoA*ez*A
        opN.M(dofU,dofA) = coupl.Tw2;  % normalized mass density
        
        % assign final radiation matrix
        opN.("R"+loading.at) = R;
        obj.opNonlin = opN;
    end
    function obj = incorporateSolidLoading(obj, loading, dofA, dofU, udof)
        % given degrees of freedom (dofA: additional for halfspace, dofU: displacements at boundary)
        opN = obj.opNonlin;
        nDof = size(opN.M,1);
        % allocate new matrices:
        Rkg = zeros(nDof); % in ik*igamma
        Rke = zeros(nDof); % in ik*ieta
        Rg  = zeros(nDof); % in igamma
        Re  = zeros(nDof); % in ieta
        % get coupling matrices: 
        coupl = PlateLeaky.couplingMatricesSolid(loading,obj.np);
        Iu = eye(3);

        % continuity of displacements ik*u - ik*ua = 0: 
        opN.L0(dofA,dofU) = -Iu(udof,udof);  % plate displacements 
        opN.L1(dofA,dofA) = +coupl.Uk(udof,udof);  % in (i*k)
        Rg(dofA,dofA) = coupl.Ug(udof,udof);      % in (i*gamma)
        Re(dofA,dofA) = coupl.Ue(udof,udof);      % in (i*eta)
        
        % balance of tractions:
        opN.L2(dofU,dofA) = opN.L2(dofU,dofA) + coupl.Tk2(udof,udof); % in (i*k)^2
        opN.M(dofU,dofA)  = opN.M(dofU,dofA)  + coupl.Tw2(udof,udof); % in w^2
        Rkg(dofU,dofA) = Rkg(dofU,dofA) + coupl.Tkg(udof,udof);    % in (i*k*i*gamma)
        Rke(dofU,dofA) = Rke(dofU,dofA) + coupl.Tke(udof,udof);    % in (i*k*i*eta)

        % assign final radiation matrix
        opN.("Rkg"+loading.at) = Rkg;
        opN.("Rke"+loading.at) = Rke;
        opN.("Rg"+loading.at) = Rg;
        opN.("Re"+loading.at) = Re;
        obj.opNonlin = opN;
    end
    function op = getPolynomialForm(obj, op, loading)
        if loading.eliminated % eliminated because beta_a = beta_b (same loading on both sides)
            return; 
        end

        side = loading.at; 
        extMat = loading.mat; 

        % solid media needs also L3 and M1: 
        if ~isfield(op,'L3') && ~hasSolidLoading(obj) 
            op.L3 = []; op.M1 = [];                % initialize if not yet done 
        elseif ~isfield(op,'L3') && hasSolidLoading(obj) 
            op.L3 = zeros(size(op.L0)); op.M1 = zeros(size(op.L0));  % initialize if not yet done 
        end

        L3 = op.L3; L2 = op.L2; L1 = op.L1; L0 = op.L0; M = op.M; M1 = op.M1;
        op = rmfield(op, {'L3','L2','L1','L0','M','M1'}); 
        Z = zeros(size(L0));

        if isa(loading.mat,"MaterialFluid")
            opName = "R"+side; 
            R = op.(opName);
            op = rmfield(op, char(opName)); 
            cf = extMat.cl/obj.np.fh0; % normalized wave speed
            Iexp = eye(2); % used to expand the remaining matrices later on 

            LL3 = kron(Iexp,L3);
            LL2 = [L2  ,   -R   ;
                   Z   ,   L2    ];
            LL1 = kron(Iexp, L1);
            LL0 = [L0  ,   Z    ; 
                   R   ,   L0   ];
            MM  = [M   ,   -1/cf^2*R  ; 
                   Z   ,      M       ];
            MM1 = kron(Iexp,M1);
        elseif isa(loading.mat,'MaterialIsotropic')
            opNames = {char("Rkg"+side), char("Rke"+side), char("Rg"+side), char("Re"+side)};
            Rkg = op.(opNames{1});
            Rke = op.(opNames{2});
            Rg  = op.(opNames{3});
            Re  = op.(opNames{4});
            op = rmfield(op, opNames); 
            cl = extMat.cl/obj.np.fh0; ct = extMat.ct/obj.np.fh0; % normalized wave speed
            al = 1/cl^2; at = 1/ct^2;
            Iexp = eye(4); % used to expand the remaining matrices later on 

            LL3 = [ L3, -Rkg,  -Rke,  Z  ; 
                    Z,   L3,     Z,  -Rke; 
                    Z,    Z,    L3,  -Rkg; 
                    Z,    Z,     Z,   L3 ];
            LL2 = [ L2, -Rg,   -Re,   Z  ; 
                    Z,   L2,    Z,    -Re; 
                    Z,    Z,   L2,    -Rg; 
                    Z,    Z,    Z,    L2]; 
            LL1 = [ L1,   Z,    Z,    Z  ; 
                    Rkg,  L1,   Z,    Z  ; 
                    Rke,  Z,    L1,   Z  ; 
                    Z,    Rke,  Rkg,  L1];
            LL0 = [ L0,   Z,    Z,    Z  ; 
                    Rg,   L0,   Z,    Z  ; 
                    Re,   Z,    L0,   Z  ; 
                    Z,    Re,   Rg,   L0];
            MM  = [  M,   -al*Rg, -at*Re,   Z  ; 
                     Z,    M,      Z,   -at*Re ;
                     Z,    Z,      M,   -al*Rg ; 
                     Z,    Z,      Z,     M   ];
            MM1 = [ M1,   -al*Rkg,   -at*Rke,   Z     ;
                    Z,       M1,        Z,    -at*Rke ;
                    Z,       Z,         M1,   -al*Rkg ; 
                    Z,       Z,         Z,      M1    ];
        else
            error('GEWTOOL:getPolynomialForm','External halfspaces need to be fluids or isotropic solids.');
        end

        % expand remaining R-matrices: 
        opList = fieldnames(op);
        for i=1:length(opList)
            Ri = op.(opList{i});
            op.(opList{i}) = kron(Iexp, Ri); % expand
        end
        op.L3 = LL3; op.L2 = LL2; op.L1 = LL1; op.L0 = LL0; op.M = MM; op.M1 = MM1;
    end
    function op = getPolynomialFormSolid(obj, op, loading)
        side = loading.at; 
        opName = "R"+side; 
        extMat = loading.mat; 
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
        % Creates a new PlateLeaky object describing only the upper symmetric half 
        % of the original object. 
        % Throws a warning if the plate is not symmetric in geometry and materials. 
        % This function is used by LambS, LambA, LambSA, etc., and you will 
        % usually not need to call it explicitly.
        if ~obj.decouplesSA('v')
            error('GEWTOOL:symmetrizeGeometry','The setup is not symmetric. I cannot symmetrize the geometry.'); 
        end
        gew = symmetrizeGeometry@PlateClosed(obj);
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
        decoupl = bothSame && decouplesSA@PlateClosed(obj,verb);
    end
    function hasSolid = hasSolidLoading(obj)
        hasSolid = false; 
        for i = 1:length(obj.halfSpaces)
            mat = obj.halfSpaces(i).mat; 
            if isa(mat,'MaterialIsotropic')
                hasSolid = true; 
                return; 
            end
        end
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
        end
        halfSpaces = struct([]); 
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
    function op = couplingMatricesSolid(loading,np)
        matA = loading.mat; % material of loading half-space
        if loading.at == "top" % different signs at top and bottom
            sig = 1; % sign
        else
            sig = -1; 
        end
        % stiffness tensor and wave velocities
        a = matA.c/np.c0; % stiffness in normalized units 
        azx = sig*squeeze(a(3,:,:,1));
        azz = sig*squeeze(a(3,:,:,3));
        cl = matA.cl/np.fh0; % longitudinal velocity in normalized units 
        ct = matA.ct/np.fh0; % transverse velocity in normalized units
        al = 1/cl^2; % w^2 ~ kappal^2 ~ 1/cl^2
        at = 1/ct^2; % w^2 ~ kappat^2 ~ 1/ct^2
        
        % coupling matrices (to be reduced to polarization "udof")
        Iu  = eye(3);
        Z   = zeros(3);
        op.Tk2  = [azx(:,1) - azz(:,3), azx(:,2), azx(:,3) + azz(:,1)]; 
        op.Tkg = [azx(:,3) + azz(:,1), Z(:,1), Z(:,1)];
        op.Tke = [Z(:,1), azz(:,2), -azx(:,1) + azz(:,3)];
        op.Tw2  = [-al*azz(:,3), Z(:,1), at*azz(:,1)]; % TODO: replace al*azz(:,3) = at*azz(:,1) = rho ??
        op.Uk  = Iu; 
        op.Ug  = [0, 0, 0 ; 0, 0, 0; 1, 0, 0]; 
        op.Ue  = [0, 0, -1; 0, 0, 0; 0, 0, 0]; 
        op.al = al; 
        op.at = at;
    end
    function op = couplingMatricesFluid(loading,np)
        matA = loading.mat; % material of loading half-space
        if loading.at == "top" % different signs at top and bottom
            sig = 1; % sign
        else
            sig = -1; 
        end
        % initialize quantities
        rhof = matA.rho/np.rho0;  % density in normalized units
        cl = matA.cl/np.fh0; % longitudinal velocity in normalized units 
        al = 1/cl^2; % w^2 ~ kappal^2 ~ 1/cl^2
        
        % coupling matrices (to be reduced to polarization "udof")
        op.Tw2 = -sig*rhof; % TODO very sign! 
        op.Ug  = 1;
        op.al = al;
    end
    function op = reduceToQuadratic(op)
        % prepare:
        if isfield(op,'L3') && isempty(op.L3)
            op = rmfield(op,'L3');
        end
        if isfield(op,'M1') && isempty(op.M1)
            op = rmfield(op,'M1');
        end
        if ~isfield(op,'L3') && ~isfield(op,'M1')
            return; 
        end

        % extract matrices and find DOFs to add: 
        L3 = op.L3; L2 = op.L2; L1 = op.L1; L0 = op.L0; M = op.M; M1 = op.M1;
        [oldRows,oldCols] = find(L3); 
        newDof = size(L3,2) + (1:length(unique(oldCols))).';  % new variables to add 
        newCols = mapOldToNew(oldCols, newDof);

        % extend all matrices:
        L3(newDof,newDof) = 0; L2(newDof,newDof) = 0; L1(newDof,newDof) = 0; L0(newDof,newDof) = 0; M(newDof,newDof) = 0; M1(newDof,newDof) = 0; 
        
        % state-space linearization:
        ind = @(r,c) sub2ind(size(L3),r,c); % find indices after expanding matrices!
        L2(ind(oldRows,newCols)) = L3(ind(oldRows,oldCols)); % linear indexing! -> not a block
         M(ind(oldRows,newCols)) = M1(ind(oldRows,oldCols)); % linear indexing! -> not a block
        L1(ind(newCols,oldCols)) = 1;  % linear indexing! -> not a block
        L0(ind(newCols,newCols)) = -1; % linear indexing! -> not a block
        % L2(ind(newCols,oldCols)) = 1;  % alternatively: reduces rank deficiency but does not regularize LL2
        % L1(ind(newCols,newCols)) = -1; % alternatively: reduces rank deficiency but does not regularize LL2
        
        % return "op":
        op = rmfield(op,{'L3','M1'});
        op.L2 = L2; op.L1 = L1; op.L0 = L0; op.M = M; 

        % local helper function:
        function newCols = mapOldToNew(oldCols, newDof)
            uc = unique(oldCols); 
            newCols = nan(size(oldCols)); 
            for i = 1:length(uc)
                newCols(oldCols == uc(i)) = newDof(i);
            end
        end
    end
end

end