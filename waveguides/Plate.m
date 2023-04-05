classdef Plate < Waveguide
% Plate - Represents guided waves in Plates.
% Displacement ansatz: u(x,y,z,t) = u(y)*exp(i k x - i w t)
% 
% Example:
% mat = Material('steel'); % load material data
% h = 1e-3; % thickness in m
% N = 20; % discretization (number of nodal points)
% plate = Plate(mat, h, N); % create waveguide description
% 
% See also Plate.Plate, Cylinder, Waveguide.
% 
% 2022 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France

methods 
	function obj = Plate(mats, ys, Ns)
        % Plate - Create a plate waveguide (propagation in x-direction)
        % Arguments: 
        % - mats:  materials [1 x Nlay], either of class "Material" or a struct
        %          needs to support mats.rho (scalar) and mats.c (3x3x3x3).
        % - ys:    coordinates of interfaces in meter [1 x Nlay+1]
        %          instead, if ys is of size [1 x Nlay]: ys are interpreted as the 
        %          thickness of each layer.
        % - Ns:    discretization order for each layer [1 x Nlay]
        % 
        % Example:
        % mat = Material('steel'); % load material data
        % h = 1e-3; % thickness in m
        % N = 20; % discretization (number of nodal points)
        % plate = Plate(mat, h, N); % create waveguide description
        % 
        if length(mats) == length(ys) % thicknesses have been provided
            if length(ys) == 1
                ys = ys*[-1/2, 1/2]; % single layer has centered coordinate 
            else
                ys = [0, cumsum(ys)]; % coordinates of interfaces starting from 0
            end
        elseif length(ys) ~= length(mats)+1
            error('GEWTOOL:Plate:wrongArguments','Provide either a thickness for each layer or the coordinates of the interfaces.');
        end
		obj = obj@Waveguide(mats, ys, Ns); % converts mats 
        obj.lay = LayerPlate.empty; % initialize with correct class
		for ii = 1:length(obj.mat)
			obj.lay(ii) = LayerPlate(obj.mat(ii), ys(ii:ii+1), Ns(ii));
		end
    end

    function gew = fullyCoupled(obj)
        n = 0;
        gew = fullyCoupled@Waveguide(obj, n);
    end

    function gew = fullyCoupledA(obj)
        % fullyCoupledA - Assemble operators for the anti-symmetric waves.
        % If your plate is not symmetric, use Plate.fullyCoupled instead.
        % 
        % See also: fullyCoupled, fullyCoupledS.
        
        obj = obj.symmetrizeGeometry;
        udof = 1:3;
		gew = obj;
		gew.op = obj.assembleLayers(udof, 0); % n = 0 (circumferential order)
        gdofs = [gew.geom.gdofBC{1}(1,1), gew.geom.gdofBC{1}(3,1)];
        gew = gew.fixGdof(gdofs); % fix ux- and uz-displacements at bottom (y=0)
    end

    function gew = fullyCoupledS(obj)
        % fullyCoupledS - Assemble operators for the symmetric waves.
        % If your plate is not symmetric, use Plate.fullyCoupled instead.
        % 
        % See also: fullyCoupled, fullyCoupledA.
        
        obj = obj.symmetrizeGeometry;
        udof = 1:3;
		gew = obj;
		gew.op = obj.assembleLayers(udof, 0); % n = 0 (circumferential order)
        gew = gew.fixGdof(gew.geom.gdofBC{1}(2,1)); % fix uy-displacement at bottom (y=0)
    end

    function gews = fullyCoupledSA(obj)
        % fullyCoupledSA - Assemble operators for the symmetric and the anti-symmetric waves.
        % The operators for symmetric (S) and anti-symmetric (A) waves are assembled separately.
        % If your plate is not symmetric, use Plate.fullyCoupled instead.
        % Return value:
        % gews: [1 x 2] array of Plate objects. 
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
    end
    
    function gew = LambS(obj)
        % LambS - Assemble operators for the symmetric Lamb polarized waves.
        % If your plate is not symmetric, use Plate.Lamb instead.
        % 
        % See also: LambA, LambSA, Lamb, sh, fullyCoupled.

        if ~obj.decouplesLambvsSH
            warning('GEWTOOL:Waveguide:donotdecouple', 'You are doing bêtises! In-plane polarized waves do not decouple from out-of plane polarization. I will proceed anyways.');
        end
        obj = obj.symmetrizeGeometry;
        udof = 1:2;
		gew = obj;
		gew.geom = Geometry(obj.geom.yItf, obj.geom.N, 2*ones(size(obj.geom.N)));
		gew.op = obj.assembleLayers(udof, 0); % n = 0 (circumferential order)
        gew = gew.fixGdof(gew.geom.gdofBC{1}(2,1)); % fix uy-displacement at bottom (y=0)
    end
    
    function gew = LambA(obj)
        % LambA - Assemble operators for the anti-symmetric Lamb polarized waves.
        % If your plate is not symmetric, use Plate.Lamb instead.
        % 
        % See also: LambS, LambSA, Lamb, sh, fullyCoupled.
        
        if ~obj.decouplesLambvsSH
            warning('GEWTOOL:Waveguide:donotdecouple', 'You are doing bêtises! In-plane polarized waves do not decouple from out-of plane polarization. I will proceed anyways.');
        end
        obj = obj.symmetrizeGeometry;
        udof = 1:2;
		gew = obj;
		gew.geom = Geometry(obj.geom.yItf, obj.geom.N, 2*ones(size(obj.geom.N)));
		gew.op = obj.assembleLayers(udof, 0); % n = 0 (circumferential order)
        gew = gew.fixGdof(gew.geom.gdofBC{1}(1,1)); % fix ux-displacement at bottom (y=0)
    end
    
    function gews = LambSA(obj)
        % LambSA - Assemble operators for the symmetric and the anti-symmetric Lamb polarized waves.
        % The operators for symmetric (S) and anti-symmetric (A) waves are assembled separately.
        % If your plate is not symmetric, use Plate.Lamb instead.
        % Return value:
        % gews: [1 x 2] array of Plate objects. 
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
            if ~m.decouplesSA
                if nargin == 2 && strcmp(verb, 'v')
                    warning('GEWTOOL:decouplesSA:matSym', 'The stiffness tensor is not invariant to reflexion ey -> -ey.');
                end
                return;
            end
        end
        % verify that the layers are distributed symmetrically:
        lMid = ceil(obj.geom.nLay/2);
        for l = 1:lMid
            if obj.lay(l) ~= obj.lay(end-(l-1))
                if nargin == 2 && strcmp(verb, 'v')
                    warning('GEWTOOL:decouplesSA:laySym', 'The layered structure is not symmetric.');
                end
                return;
            end
        end
        decoupl = true;
    end
    
    function gew = symmetrizeGeometry(obj)
        % symmetrizeGeometry - upper symmetric half of the original plate.
        % Creates a new Plate object describing only the upper symmetric half 
        % of the original Plate object. 
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
        yItfList = [obj.geom.yItf(lMid:end,1).' obj.geom.yItf(end,end)];
        if mod(obj.geom.nLay,2) == 1 % odd number of layers
            hmid = yItfList(2) - yItfList(1);
            yItfList(1) = yItfList(2) - hmid/2; % half the thickness for middle layer
        end
        yItfList = yItfList - yItfList(1); % zero coordinate at center
        gew = Plate([lays.mat], yItfList, Ns);
    end

end % methods

end % class
