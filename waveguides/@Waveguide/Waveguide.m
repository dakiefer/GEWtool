classdef Waveguide < matlab.mixin.Copyable

properties (Access = public)
	geom
	mat
	lay Layer = LayerPlate.empty
	op = [] 	% operators 
	np  		% normalization parameters (material, geom)
end % properties

properties (Dependent)
	h0  		% short hand for obj.np.h0
end % properties Dependent


methods 
	function obj = Waveguide(mats, rs, Ns, Nudof)
		if nargin < 4, Nudof = 3*ones(size(Ns)); end
		obj.geom = Geometry(rs, Ns, Nudof);
		obj.mat = mats;
		np.c0 = mats(1).c(1,2,1,2); 
		np.h0 =(rs(end)-rs(1))/length(mats); % normalization parameters
		np.rho0 = mats(1).rho;
		np.fh0 = sqrt(np.c0/np.rho0);
		obj.np = np;
	end

	function h0 = get.h0(obj)
		h0 = obj.np.h0;
	end

	function guw = fullyCoupled(obj, n)
		udof = 1:3;
		guw = obj;
        guw.geom = Geometry(obj.geom.yItf, obj.geom.N, 3*ones(size(obj.geom.N))); % update 
		guw.op = obj.assembleLayers(udof, n);
		guw.op = obj.freeBCs(udof, n);
	end

	function guw = Lamb(obj, n)
		udof = 1:2;
		guw = obj;
		guw.geom = Geometry(obj.geom.yItf, obj.geom.N, 2*ones(size(obj.geom.N)));
		guw.op = obj.assembleLayers(udof, n);
		guw.op = obj.freeBCs(udof, n);
	end

	function obj = sh(obj, n)
		udof = 3;
		guw = obj;
		guw.geom = Geometry(obj.geom.yItf, obj.geom.N, ones(size(obj.geom.N)));
		guw.op = obj.assembleLayers(udof, n);
		guw.op = obj.freeBCs(udof, n);
	end

	[op] = assembleLayers(obj, udof, n)
	[op] = freeBCs(obj, udof, n)

end

end % class
