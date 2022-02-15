classdef Plate < Waveguide

methods 
	function obj = Plate(mats, ys, Ns)
        if isscalar(ys), ys = ys*[-1/2, 1/2]; end % convert thickness to y-range
		obj = obj@Waveguide(mats, ys, Ns);
		for ii = 1:length(mats)
			obj.lay(ii) = LayerPlate(mats(ii), ys(ii:ii+1), Ns(ii));
		end
	end
end % methods

end % class
