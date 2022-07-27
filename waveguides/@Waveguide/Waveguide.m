classdef Waveguide < matlab.mixin.Copyable
% Waveguide - Represents guided waves in plates or cylinders.
% Usually, you do not need to use this class directly. For simpler interfacing, use 
% the derived classes "Plate" or "Cylinder" instead.
% 
% See also Plate, Cylinder.
% 
% 2022 - Daniel A. Kiefer
% Institut Langevin, Paris, France
% 

properties (Access = public)
	geom        % geometry object describing the discretized, multilayered structure
	mat         % materials [1 x Nlay]
	lay Layer = LayerPlate.empty   % layers [1 x Nlay]
	op = [] 	% operators describing the wave propagation
	np  		% normalization parameters (material parameters, geometry)
end % properties

properties (Dependent)
	h0  		% unit thickness: short hand for obj.np.h0
end % properties Dependent


methods 
	function obj = Waveguide(mats, rs, Ns, Nudof)
        % Waveguide - Create a waveguide object.
        % Arguments: 
        % - mats:  materials [1 x Nlay], either of class "Material" or a struct
        %          needs to support mats.rho (scalar) and mats.c (3x3x3x3).
        % - rs:    coordinates of interfaces in meter [1 x Nlay+1]
        % - Ns:    discretization order [1 x Nlay], each entry corresponds to one layer
        % - Nudof: displacement digrees of freedom [1 x Nlay], each entry corresponds to one layer
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

	function gew = fullyCoupled(obj, n)
        % fullyCoupled - Assemble wave operators describing the coupled set of Lamb- and SH-polarized waves.
        % 
        % See also: Lamb, sh, decouplesLambvsSH.
		udof = 1:3;
		gew = obj;
        gew.geom = Geometry(obj.geom.yItf, obj.geom.N, 3*ones(size(obj.geom.N))); % update 
		gew.op = obj.assembleLayers(udof, n);
	end

	function gew = Lamb(obj, n)
        % Lamb - Assemble wave operators describing the Lamb polarized waves (in-plane).
        % 
        % See also: sh, fullyCoupled, decouplesLambvsSH.
        if ~obj.decouplesLambvsSH
            warning('GEWTOOL:Waveguide:donotdecouple', 'You are doing bêtises! In-plane polarized waves do not decouple from out-of plane polarization. I will proceed anyways.');
        end
		udof = 1:2;
		gew = obj;
		gew.geom = Geometry(obj.geom.yItf, obj.geom.N, 2*ones(size(obj.geom.N)));
		gew.op = obj.assembleLayers(udof, n);
    end

	function gew = sh(obj, n)
        % sh - Assemble wave operators describing the shear-horizontal (SH) polarized waves (out-of-plane).
        %
        % See also: Lamb, fullyCoupled, decouplesLambvsSH.
        if ~obj.decouplesLambvsSH
            warning('GEWTOOL:Waveguide:donotdecouple', 'You are doing bêtises! In-plane polarized waves do not decouple from out-of plane polarization. I will proceed anyways.');
        end
		udof = 3;
		gew = obj;
		gew.geom = Geometry(obj.geom.yItf, obj.geom.N, ones(size(obj.geom.N)));
		gew.op = obj.assembleLayers(udof, n);
    end
    
    function lin = isLinearizableInK2(obj)
        % isLinearizableInK2 - Test if the problem can be linearized in k without increasing the problem size.
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
        % linearizeInK2 - Linearizes the problem in the wavenumber k without increasing size. 
        % The matrices are manipulated without increasing their size. Use this before passing 
        % to the solver for faster computation of wavenumbers.
        if ~obj.isLinearizableInK2
            warning('GEWTOOL:Waveguide:nonlinearizable', 'This problem is not linearizable as intended.');
        end
        L2 = obj.op.L2; L1 = obj.op.L1; L0 = obj.op.L0;
        N = obj.geom.N;
        dofx = 1:N; dofy = N+1:2*N;
        L2(dofy,dofx) = L1(dofy,dofx); 
        L0(dofx,dofy) = L1(dofx,dofy);
        obj.op.L2 = L2; obj.op.L1 = []; obj.op.L0 = L0;
    end
    
    function decoupl = decouplesLambvsSH(obj)
        % decouplesLambvsSH - Tests whether the Lamb- and SH-polarized waves decouple.
        % 
        % See also: Lamb, sh, fullyCoupled.
        lays = obj.lay;
        for l = lays % test each of the layers
            if ~l.decouplesLambvsSH
                decoupl = false; return;
            end
        end
        decoupl = true;
    end
    
    function obj = fixGdof(obj, gdof)
        % fixGdof - Homogeneous Dirichlet BCs: Fixes the specified degrees of freedom.
        % The displacements at the indicated global degrees of freedom are set to zero.
        % This function is used internally to implement the symmetric and anti-symmetric
        % Lamb waves.
        if isempty(obj.op)
            warning('GEWTOOL:Waveguide:notassembled', 'Define the waveguide problem first by calling, e.g., fullyCoupled().');
            return
        end
        ops = fields(obj.op);
        for i=1:length(ops)
            opName = ops{i};
            obj.op.(opName)(gdof,:) = []; % remove row
            obj.op.(opName)(:,gdof) = []; % remove column
        end
        for l=1:obj.geom.nLay
%             obj.geom.gdofOfLay{l} = setdiff(obj.geom.gdofOfLay{l}, gdof);
%             obj.geom.gdofBC{l} = setdiff(obj.geom.gdofBC{l}, gdof);
        end
        obj.geom.gdofDBC = [obj.geom.gdofDBC(:) gdof(:).'];
    end

	[op] = assembleLayers(obj, udof, n)
end

end % class
