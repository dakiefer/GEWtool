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

    function guw = sh(obj)
        n = 0;
        guw = sh@Waveguide(obj, n);
    end

end % methods

end % class
