classdef PlateClosed < Waveguide
% PlateClosed - Represents guided waves in Plates.
% Displacement ansatz: u(x,y,z,t) = u(z)*exp(i k x - i w t)
% 
% Example:
% mat = Material('steel'); % load material data
% h = 1e-3; % thickness in m
% N = 20; % discretization (number of nodal points)
% plate = PlateClosed(mat, h, N); % create waveguide description
% 
% See also PlateClosed.PlateClosed, Cylinder, Waveguide.
% 
% 2022 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France

methods 
	function obj = PlateClosed(mats, zs, Ns)
        % PlateClosed - Create a plate waveguide (propagation in x-direction)
        % Arguments: 
        % - mats:  materials [1 x Nlay], either of class "Material" or a struct
        %          needs to support mats.rho (scalar) and mats.c (3x3x3x3).
        % - zs:    coordinates of interfaces in meter [1 x Nlay+1]
        %          instead, if zs is of size [1 x Nlay]: zs are interpreted as the 
        %          thickness of each layer.
        % - Ns:    discretization order for each layer [1 x Nlay]
        % 
        % Example:
        % mat = Material('steel'); % load material data
        % h = 1e-3; % thickness in m
        % N = 20; % discretization (number of nodal points)
        % plate = PlateClosed(mat, h, N); % create waveguide description
        % 
		obj = obj@Waveguide(mats, zs, Ns); % converts mats 
		for ii = 1:length(obj.mat)
            if isa(obj.mat{ii}, 'MaterialPiezoelectric') % subclass first! 
                obj.lay{ii} = LayerPlatePiezo(obj.mat{ii}, obj.geom.zItf(ii,:), obj.geom.N(ii));
            elseif isa(obj.mat{ii}, 'Material') % subclasses are also Material objects
			    obj.lay{ii} = LayerPlate(obj.mat{ii}, obj.geom.zItf(ii,:), obj.geom.N(ii));
            else 
                error('GEWtool:PlateClosed', ['The material type "%s" is not ' ...
                    'known. It should be of class "Material", "MaterialPiezoelectric" ' ...
                    'or a subclass thereof.'],class(obj.mat{ii}) );
            end
		end
    end

    function gew = fullyCoupled(obj)
        n = 0;
        gew = fullyCoupled@Waveguide(obj, n);
        [gew.family] = deal('all'); % assign same to each element of gew if necessary
    end

    function gew = fullyCoupledA(obj)
        % fullyCoupledA - Assemble operators for the anti-symmetric waves.
        % If your plate is not symmetric, use PlateClosed.fullyCoupled instead.
        % 
        % See also: fullyCoupled, fullyCoupledS.
        
        if ~obj.geom.symmetrized, obj = obj.symmetrizeGeometry; end
        udof = 1:3;
        udofFix = obj.dofInplane(udof); 
        gew = obj.polarization(udof,0); % n = 0 (circumferential order)
        gew = gew.addDOFtoDBC(gew.geom.gdofBC{1}(udofFix,1)); % fix ux- and uy-displacements at bottom (z=0)
        gew = gew.assembleLayers(udof,0);
        [gew.family] = deal('all anti-sym.'); % assign same to each element of gew if necessary
    end

    function gew = fullyCoupledS(obj)
        % fullyCoupledS - Assemble operators for the symmetric waves.
        % If your plate is not symmetric, use PlateClosed.fullyCoupled instead.
        % 
        % See also: fullyCoupled, fullyCoupledA.
        
        if ~obj.geom.symmetrized, obj = obj.symmetrizeGeometry; end
        udof = 1:3;
        udofFix = obj.dofOutofplane(udof); % 3 if [ux, uy, uz], 2 if [ux, uz]
		gew = obj.polarization(udof,0); % n = 0 (circumferential order)
        gew = gew.addDOFtoDBC(gew.geom.gdofBC{1}(udofFix,1)); % fix uz-displacement at bottom (z=0)
        gew = gew.assembleLayers(udof,0);
        [gew.family] = deal('all sym.'); % assign same to each element of gew if necessary
    end

    function gews = fullyCoupledSA(obj)
        % fullyCoupledSA - Assemble operators for the symmetric and the anti-symmetric waves.
        % The operators for symmetric (S) and anti-symmetric (A) waves are assembled separately.
        % If your plate is not symmetric, use PlateClosed.fullyCoupled instead.
        % Return value:
        % gews: [1 x 2] array of PlateClosed objects. 
        %       - gews(1) describes the symmetric waves
        %       - gews(2) describes the anti-symmetric waves
        % 
        % See also: fullyCoupled, fullyCoupledA, fullyCoupledS, Lamb, sh.

        gews(1) = obj.fullyCoupledS;
        gews(2) = obj.fullyCoupledA;
    end

    function gew = Lamb(obj)
        n = 0;
        gew = Lamb@Waveguide(obj, n);
        [gew.family] = deal('Lamb'); % assign same to each element of gew if necessary
    end
    
    function gew = LambS(obj)
        % LambS - Assemble operators for the symmetric Lamb polarized waves.
        % If your plate is not symmetric, use PlateClosed.Lamb instead.
        % 
        % See also: LambA, LambSA, Lamb, sh, fullyCoupled.

        if ~obj.geom.symmetrized, obj = obj.symmetrizeGeometry; end
        udof = obj.udofLamb;
        udofFix = obj.dofOutofplane(udof); % 3 if [ux, uy, uz], 2 if [ux, uz]
		gew = obj.polarization(udof,0); % n = 0 (circumferential order)
        gew = gew.addDOFtoDBC(gew.geom.gdofBC{1}(udofFix,1)); % fix uz-displacement at bottom (z=0)
        gew = gew.assembleLayers(udof,0);
        [gew.family] = deal('Lamb sym.'); % assign same to each element of gew if necessary
    end
    
    function gew = LambA(obj)
        % LambA - Assemble operators for the anti-symmetric Lamb polarized waves.
        % If your plate is not symmetric, use PlateClosed.Lamb instead.
        % 
        % See also: LambS, LambSA, Lamb, sh, fullyCoupled.
        
        if ~obj.geom.symmetrized, obj = obj.symmetrizeGeometry; end
        udof = obj.udofLamb;
        udofFix = obj.dofInplane(udof);
		gew = obj.polarization(udof,0); % n = 0 (circumferential order)
        gew = gew.addDOFtoDBC(gew.geom.gdofBC{1}(udofFix,1)); % fix ux-displacement at bottom (z=0)
        gew = gew.assembleLayers(udof,0);
        [gew.family] = deal('Lamb anti-sym.'); % assign same to each element of gew if necessary
    end
    
    function gews = LambSA(obj)
        % LambSA - Assemble operators for the symmetric and the anti-symmetric Lamb polarized waves.
        % The operators for symmetric (S) and anti-symmetric (A) waves are assembled separately.
        % If your plate is not symmetric, use PlateClosed.Lamb instead.
        % Return value:
        % gews: [1 x 2] array of PlateClosed objects. 
        %       - gews(1) describes the symmetric waves
        %       - gews(2) describes the anti-symmetric waves
        % 
        % See also: Lamb, LambA, LambS, sh, fullyCoupled.

        gews(1) = obj.LambS;
        gews(2) = obj.LambA;
    end

    function gew = sh(obj)
        n = 0;
        gew = sh@Waveguide(obj, n);
        [gew.family] = deal('SH waves'); % assign same to each element of gew if necessary
    end

    function decoupl = decouplesSA(obj, verb)
        % decouplesSA - Tests whether symmetric and antisymmetric waves decouple.
        % Usage: 
        % decoupl = decouplesSA;       Returns true if SA waves decouple.
        % decoupl = decouplesSA('v');  Throw warning indicating reason.
        % 
        % See also: decouplesLambvsSH.

        % verify symmetry of materials:
        mats = obj.mat; % array of materials for each layer
        decoupl = false; 
        for m = mats
            if ~m{1}.decouplesSA
                if nargin == 2 && strcmp(verb, 'v')
                    warning('GEWTOOL:decouplesSA:matSym', 'The stiffness tensor is not invariant to reflexion ez -> -ez.');
                end
                return;
            end
        end
        % verify that the layers are distributed symmetrically:
        lMid = ceil(obj.geom.nLay/2);
        for l = 1:lMid
            if obj.lay{l} ~= obj.lay{end-(l-1)}
                error('GEWTOOL:decouplesSA:laySym', 'The layered structure is not symmetric.');
            end
        end
        decoupl = true;
    end
    
    function gew = symmetrizeGeometry(obj)
        % symmetrizeGeometry - upper symmetric half of the original plate.
        % Creates a new PlateClosed object describing only the upper symmetric half 
        % of the original PlateClosed object. 
        % Throws a warning if the plate is not symmetric in geometry and materials. 
        % This function is used by LambS, LambA, LambSA, etc., and you will 
        % usually not need to call it explicitly.

        if ~obj.decouplesSA('v') % 'v' verbose: throw warning indicating reason
            warning('GEWTOOL:symmetrizeGeometry:SAdecoupl', 'You are doing bêtises! S/A waves do not decouple. I will proceed anyways.');
        end
        % create new waveguide object reduced to only the symmetric half:
        if mod(obj.geom.nLay,2) == 1 % odd number of layers
            lMid = ceil(obj.geom.nLay/2);
        else % even number of layers
            lMid = obj.geom.nLay/2+1;
        end
        lays = obj.lay(lMid:end);
        Ns = obj.geom.N(lMid:end);
        Ns = ceil(Ns/2); % use smaller matrices
        while any(Ns<2), Ns = Ns +1; end % ensure to have at least two nodes
        zItfList = [obj.geom.zItf(lMid:end,1).' obj.geom.zItf(end,end)];
        if mod(obj.geom.nLay,2) == 1 % odd number of layers
            hmid = zItfList(2) - zItfList(1);
            zItfList(1) = zItfList(2) - hmid/2; % half the thickness for middle layer
        end
        zItfList = zItfList - zItfList(1); % zero coordinate at center
        mats = cell(1,length(lays)); % allocate
        for i = 1:length(lays), mats{i} = lays{i}.mat; end % extract materials 
        constructorFun = str2func(class(obj)); % instantiate of same class as obj
        gew = constructorFun(mats, zItfList, Ns);
        gew.geom.symmetrized = true;
    end

    function F = displGrad(obj, dat)
        % displGrad - Displacement gradient of the provided field.
        F = cell(obj.geom.nLay, 1); % allocate for each layer
        u = displacement(dat); % extract displacements (piezoelec. plate also has potentials)
        for l = 1:obj.geom.nLay
            NudofL = size(u{l},4);
            sizeF = [size(dat.w), obj.geom.N(l), NudofL, NudofL]; % size of F = grad u
            Fi = zeros(sizeF); % allocate for displacement gradient F = grad u
            lay = obj.lay{l};
            Dz = 1/lay.h*lay.D1*obj.np.h0;  % differentiation matrix, normalized
            iku = 1i*dat.k.*u{l}*obj.np.h0; % ex.F = ik*u, normalized
            uu = permute(u{l}, [1, 2, 5, 3, 4]); % additional dimension for mult. with diff. mat.
            dzu = sum(shiftdim(Dz, -2).*uu, 4); % ez.F = ∂u/∂z, dimension 4 is singleton
            Fi(:,:,:,1,:) = iku;  % assign components ex.F
            Fi(:,:,:,NudofL,:) = dzu;  % assign components ez.F
            F{l} = Fi; 
        end 
    end

    function G = potentialGrad(obj, dat)
        % potentialGrad - Gradient of electric potential provided in "dat"
        % G = grad(phi) = ex ik phi + ez ∂z phi
        G = cell(obj.geom.nLay, 1); % allocate for each layer
        phi = potential(dat);       % extract potentials (zero if no piezoelectricity)
        Nudof = length(obj.udof);   % number of displacement components
        xDim = 1;                  % index of x-components
        zDim = obj.dofOutofplane(obj.udof); % index of z-components (2 or 3, depending on polarization)
        for l = 1:obj.geom.nLay
            sizeG = [size(dat.w), obj.geom.N(l), Nudof]; % same size as displacement vectors (depends on polarization)
            Gi = zeros(sizeG); % allocate for this layer
            if isPiezoelectric(dat.gew)
                lay = obj.lay{l};
                Dz = 1/lay.h*lay.D1*obj.np.h0; % differentiation matrix, normalized
                pphi = permute(phi{l}, [1, 2, 4, 3]); % additional dimension for mult. with diff. mat.
                dzPhi = sum(shiftdim(Dz, -2).*pphi, 4); % ez.F = ∂u/∂z, dimension 4 is singleton
                Gi(:,:,:,xDim,:) = 1i*(dat.k*obj.np.h0).*phi{l};  % assign components ex.F, normalized
                Gi(:,:,:,zDim,:) = dzPhi;  % assign components ez.F, already normalized
            end
            G{l} = Gi; 
        end 
    end

end % methods

end % class
