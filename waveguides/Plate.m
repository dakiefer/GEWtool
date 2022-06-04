classdef Plate < Waveguide

methods 
	function obj = Plate(mats, ys, Ns)
        if length(mats) == length(ys) % thicknesses have been provided
            if length(ys) == 1
                ys = ys*[-1/2, 1/2]; % single layer has centered coordinate 
            else
                ys = [0, cumsum(ys)]; % coordinates of interfaces starting from 0
            end
        elseif length(ys) ~= length(mats)+1
            error('Provide either a thickness for each layer or the coordinates of the interfaces.');
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
        n = 0;
        obj = obj.symmetrizeGeometry;
        guw = LambS@Waveguide(obj, n);
    end
    
    function guw = LambA(obj)
        n = 0;
        obj = obj.symmetrizeGeometry;
        guw = LambA@Waveguide(obj, n);
    end
    
    function guws = LambSA(obj)
        guws(1) = obj.LambS;
        guws(2) = obj.LambA;
    end

    function guw = sh(obj)
        n = 0;
        guw = sh@Waveguide(obj, n);
    end
    
    function guw = symmetrizeGeometry(obj)
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
