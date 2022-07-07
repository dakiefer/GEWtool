classdef Plate < Waveguide
% Plate - Represents guided waves in Plates.
% Mechanical displacement ansatz: u(x,y,z,t) = u(y)*exp(i k x - i w t)
% 
% See also Plate.Plate, Cylinder, Waveguide.
% 
% 2022 - Daniel A. Kiefer
% Institut Langevin, Paris, France
% 

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
		obj = obj@Waveguide(mats, ys, Ns);
        obj.lay = LayerPlate.empty; % initialize with correct class
		for ii = 1:length(mats)
			obj.lay(ii) = LayerPlate(mats(ii), ys(ii:ii+1), Ns(ii));
		end
    end

    function guw = fullyCoupled(obj)
        n = 0;
        guw = fullyCoupled@Waveguide(obj, n);
    end

    function guw = Lamb(obj)
        n = 0;
        guw = Lamb@Waveguide(obj, n);
    end
    
    function guw = LambS(obj)
        % LambS - Assemble operators for the symmetric Lamb polarized waves.
        % If your plate is not symmetric, use Plate.Lamb instead.
        % 
        % See also: LambA, LambSA, Lamb, sh, fullyCoupled.
        n = 0;
        obj = obj.symmetrizeGeometry;
        udof = 1:2;
		guw = obj;
		guw.geom = Geometry(obj.geom.yItf, obj.geom.N, 2*ones(size(obj.geom.N)));
		guw.op = obj.assembleLayers(udof, n);
        guw = guw.fixGdof(guw.geom.gdofBC{1}(2,1)); % fix uy-displacement at bottom (y=0)
    end
    
    function guw = LambA(obj)
        % LambA - Assemble operators for the anti-symmetric Lamb polarized waves.
        % If your plate is not symmetric, use Plate.Lamb instead.
        % 
        % See also: LambS, LambSA, Lamb, sh, fullyCoupled.
        n = 0;
        obj = obj.symmetrizeGeometry;
        udof = 1:2;
		guw = obj;
		guw.geom = Geometry(obj.geom.yItf, obj.geom.N, 2*ones(size(obj.geom.N)));
		guw.op = obj.assembleLayers(udof, n);
        guw = guw.fixGdof(guw.geom.gdofBC{1}(1,1)); % fix ux-displacement at bottom (y=0)
    end
    
    function guws = LambSA(obj)
        % LambSA - Assemble operators for the symmetric and the anti-symmetric Lamb polarized waves.
        % The operators for symmetric (S) and anti-symmetric (A) waves are assembled separately.
        % If your plate is not symmetric, use Plate.Lamb instead.
        % Return value:
        % guws: [1 x 2] array of Plate objects. 
        %       - guws(1) describes the anti-symmetric waves
        %       - guws(2) describes the symmetric waves
        % 
        % See also: Lamb, LambA, LambS, sh, fullyCoupled.
        guws(1) = obj.LambS;
        guws(2) = obj.LambA;
    end

    function guw = sh(obj)
        n = 0;
        guw = sh@Waveguide(obj, n);
    end
    
    function guw = symmetrizeGeometry(obj)
        % symmetrizeGeometry - upper symmetric half of the original plate.
        % Creates a new Plate object describing only the upper symmetric half 
        % of the original Plate object. The input plate should be symmetric in 
        % geometry and materials. This function is used by LambS, LambA, LambSA
        % and you will usually not need to call it explicitly.

        % TODO verify symmetry
        lMid = ceil(obj.geom.nLay/2);
        lays = obj.lay(lMid:end);
        Ns = obj.geom.N(lMid:end);
        Ns = ceil(Ns/2); % use smaller matrices
        yItfList = [obj.geom.yItf(lMid:end,1).' obj.geom.yItf(end,end)];
        if mod(obj.geom.nLay,2) == 1
            hmid = yItfList(2) - yItfList(1);
            yItfList(1) = yItfList(2) - hmid/2; % half the thickness for middle layer
        end
        yItfList = yItfList - yItfList(1); % zero coordinate at center
        guw = Plate([lays.mat], yItfList, Ns);
    end

end % methods

end % class
