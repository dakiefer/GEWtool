classdef Cylinder < Waveguide

methods 
	function obj = Cylinder(mats, rs, Ns)
		obj = obj@Waveguide(mats, rs, Ns);
        obj.lay = LayerCylindrical.empty; % initialize with correct class
		for ii = 1:length(mats)
			obj.lay(ii) = LayerCylindrical(mats(ii), rs(ii:ii+1), Ns(ii));
		end
	end
end % methods

end % class