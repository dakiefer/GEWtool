classdef Plate < Waveguide

methods 
	function obj = Plate(mats, ys, Ns)
        if isscalar(ys) 
            ys = ys*[-1/2, 1/2]; % convert thickness to y-range
        elseif length(ys) == length(mats)
            ys = [0; cumsum(ys(:))].';
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
