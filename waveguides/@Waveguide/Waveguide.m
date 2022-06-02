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
	end

	function guw = Lamb(obj, n)
        if ~obj.decouplesLambvsSH
            warning('GUWtool:Waveguide You are doing bêtises! In-plane polarized waves do not decouple from out-of plane polarization. I will proceed anyways.');
        end
		udof = 1:2;
		guw = obj;
		guw.geom = Geometry(obj.geom.yItf, obj.geom.N, 2*ones(size(obj.geom.N)));
		guw.op = obj.assembleLayers(udof, n);
	end

	function obj = sh(obj, n)
        if ~obj.decouplesLambvsSH
            warning('GUWtool:Waveguide You are doing bêtises! In-plane polarized waves do not decouple from out-of plane polarization. I will proceed anyways.');
        end
		udof = 3;
		guw = obj;
		guw.geom = Geometry(obj.geom.yItf, obj.geom.N, ones(size(obj.geom.N)));
		guw.op = obj.assembleLayers(udof, n);
    end
    
    function lin = isLinearizableInK2(obj)
        if obj.geom.Nudof < 2, lin = false; return; end % TODO look at SH waves
        if obj.geom.nLay > 1, lin = false; return; end  % TODO implement for multiple layers
        if isempty(obj.op), warning('Setup problem first.'); end
        N = obj.geom.N; 
        dofx = 1:N; dofy = N+1:2*N;
        L2test = all(obj.op.L2(dofy,dofx) == 0, 'all');
        L1test = all(obj.op.L1(dofx,dofx) == 0, 'all') & all(obj.op.L1(dofy,dofy) == 0, 'all');
        L0test = all(obj.op.L0(dofx,dofy) == 0, 'all');
        lin = L2test & L1test & L0test;
    end
    
    function obj = linearizeInK2(obj)
        L2 = obj.op.L2; L1 = obj.op.L1; L0 = obj.op.L0;
        N = obj.geom.N;
        dofx = 1:N; dofy = N+1:2*N;
        L2(dofy,dofx) = L1(dofy,dofx); 
        L0(dofx,dofy) = L1(dofx,dofy);
        obj.op.L2 = L2; obj.op.L1 = []; obj.op.L0 = L0;
    end
    
    function decoupl = decouplesLambvsSH(obj)
        lays = obj.lay;
        for l = lays % test each of the layers
            if ~l.decouplesLambvsSH
                decoupl = false; return;
            end
        end
        decoupl = true;
    end

	[op] = assembleLayers(obj, udof, n)
	[op] = freeBCs(obj, udof, n)

end

end % class
