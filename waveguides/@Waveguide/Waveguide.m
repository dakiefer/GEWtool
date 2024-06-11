classdef Waveguide < matlab.mixin.Copyable
% Waveguide - Represents guided waves in plates or cylinders.
% Waveguide objects are handle objects, i.e., they are passed by reference.
% Usually, you do not need to use this class directly. For simpler interfacing, use 
% the derived classes "Plate" or "Cylinder" instead.
% 
% See also Plate, Cylinder.
% 
% 2022 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France

properties (Access = public)
	geom        % geometry object describing the discretized, multilayered structure
	lay Layer = LayerPlate.empty   % layers [1 x Nlay]
	op = [] 	% operators describing the wave propagation
    udof = []   % polarization: displacement components accounted for
	np  		% normalization parameters (material parameters, geometry)
end % properties

properties (Access = protected)
    mat         % materials [1 x Nlay] (user interface should rather use obj.lay.mat)
end

properties (Dependent)
    h           % total thickness
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

        % initialize geometry:
        if isscalar(Ns) && length(mats) > 1 % use same discretization for all layers
            Ns = Ns*ones(size(mats)); % expand into a vector of nodal numbers
        end
		if nargin < 4, Nudof = 3*ones(size(Ns)); end % do after expanding Ns!
        if any(diff(rs) <= 0)
            error('GEWTOOL:Waveguide',...
            'All layers should have a positive thickness.');
        end
        obj.geom = Geometry(rs, Ns, Nudof);
        % convert mats to Material class:
        if isstruct(mats) % convert from struct to objects of Material class
            matsObj = Material.empty(0, length(mats)); % initialize
            for i = 1:length(mats)
                matsObj(i) = Material(mats(i));
            end
            mats = matsObj;
        end
		obj.mat = mats; % protected property is later used in constructor of subclass
        % choose normalization parameters (physical units for the calculation):
        c1212 = zeros(1,length(mats));
        for i = 1:length(mats), c1212(i) = mats(i).c(1,2,1,2); end
		np.c0 = mean(c1212); % unit stiffness
		np.rho0 = mean([mats.rho]); % unit mass
        np.h0 = obj.h/length(mats); % unit distance
		np.fh0 = sqrt(np.c0/np.rho0); % unit frequency-thickness
		obj.np = np; % normalization parameters
	end

    function h = get.h(obj)
        h = obj.geom.yItf(end)-obj.geom.yItf(1);
    end

    function obj = polarization(obj, udof, n)
        % POLARIZATION - Assemble wave operators for given polarization and order.
        % Argument:
        % - udof:   desired polarization (displacement components) as a vector.
        %           e.g. [1 2 3], [1 2], [3].
        % - n:      order of the waves (circumferential order in cylinders).
        %           ignore if not needed. default: 0.
        % 
        % See also: Lamb, sh, decouplesLambvsSH.
        Nudof = length(udof);
        if any(obj.geom.Nudof ~= Nudof) % update geometry if necessary
            geomNew = Geometry(obj.geom.yItf, obj.geom.N, Nudof*ones(obj.geom.nLay,1)); 
            geomNew.symmetrized = obj.geom.symmetrized;
            obj.geom = geomNew; 
        end
		obj.assembleLayers(udof, n);
        obj.udof = udof;  % remember polarization
        if isa(obj,'Cylinder'), obj.n = n; end  % remember circumferential wavenumber
	end

	function gew = fullyCoupled(obj, n)
        % fullyCoupled - Assemble wave operators describing the coupled set of Lamb- and SH-polarized waves.
        % 
        % See also: Lamb, sh, decouplesLambvsSH.
		udof = 1:3;
        gew = obj.polarization(udof, n);
	end

	function gew = Lamb(obj, n)
        % Lamb - Assemble wave operators describing the Lamb polarized waves (in-plane).
        % 
        % See also: sh, fullyCoupled, decouplesLambvsSH.
        if ~obj.decouplesLambvsSH(n)
            warning('GEWTOOL:Waveguide:donotdecouple', 'You are doing bêtises! In-plane polarized waves do not decouple from out-of plane polarization. I will proceed anyways.');
        end
		udof = 1:2;
        gew = obj.polarization(udof, n);
    end

	function gew = sh(obj, n)
        % sh - Assemble wave operators describing the shear-horizontal (SH) polarized waves (out-of-plane).
        %
        % See also: Lamb, fullyCoupled, decouplesLambvsSH.
        if ~obj.decouplesLambvsSH(n)
            warning('GEWTOOL:Waveguide:donotdecouple', 'You are doing bêtises! In-plane polarized waves do not decouple from out-of plane polarization. I will proceed anyways.');
        end
		udof = 3;
		gew = obj.polarization(udof, n);
    end

    function dis = isDissipative(obj)
        % ISDISSIPATIVE - Check if waveguide dissipates energy. 
        % returns
        % - true: if any of the layers is composed of a dissipative material
        % - false: otherwise
        dis = true; 
        for l = obj.lay
            if ~isDissipative(l.mat)
                dis = false; 
                return;
            end
        end
    end
    
    function lin = isLinearizableInK(obj)
        % isLinearizableInK - Test if the problem can be linearized in k without increasing the problem size.
        if isempty(obj.op)
            warning('GEWTOOL:isLinearizableInK', 'Setup problem first. Ignoring your request.'); 
            lin = false; return; 
        end
        if isempty(obj.op.L1) && ~isempty(obj.op.L2)
            lin = true; return; % already linearized
        end
        if obj.geom.Nudof < 2, lin = false; return; end % TODO look at SH waves
        if obj.geom.nLay > 1, lin = false; warning('Not implemented for multilayers.'); return; end  % TODO implement for multiple layers
        dofx = obj.geom.gdofRedX; dofy = obj.geom.gdofRedY;
        L2test = all(obj.op.L2(dofy,dofx) == 0, 'all');
        L1test = all(obj.op.L1(dofx,dofx) == 0, 'all') & all(obj.op.L1(dofy,dofy) == 0, 'all');
        L0test = all(obj.op.L0(dofx,dofy) == 0, 'all');
        lin = L2test & L1test & L0test;
    end
    
    function gews = linearizeInK(gews)
        % linearizeInK - Linearizes the problem in the wavenumber k without increasing size. 
        % The matrices are manipulated without increasing their size. Use this before passing 
        % to the solver for faster computation of wavenumbers.
        %
        % Literature
        % [1] E. Kausel, “An explicit solution for the Green functions for
        % dynamic loads in layered media,” Massachusetts Institute Of
        % Technology, Research Report R81-13, 1981.
        % [2] H. Gravenkamp, “Efficient simulation of elastic guided waves
        % interacting with notches, adhesive joints, delaminations and inclined
        % edges in plate structures,” Ultrasonics, vol. 82, pp. 101–113, Jan.
        % 2018, doi: 10.1016/j.ultras.2017.07.019.
        for obj = gews % note: the objects are "by reference", i.e., the original ones are changed
            if ~obj.isLinearizableInK
                warning('GEWTOOL:Waveguide:nonlinearizable', 'This problem is not linearizable as intended.');
            end
            if ~isempty(obj.op.L1) % ignore if already linearized
                L2 = obj.op.L2; L1 = obj.op.L1; L0 = obj.op.L0;
                dofx = obj.geom.gdofRedX; dofy = obj.geom.gdofRedY;
                L2(dofy,dofx) = L1(dofy,dofx); 
                L0(dofx,dofy) = L1(dofx,dofy);
                obj.op.L2 = L2; obj.op.L1 = []; obj.op.L0 = L0;
            end
        end
    end
    
    function decoupl = decouplesLambvsSH(obj,n)
        % decouplesLambvsSH - Tests whether the Lamb- and SH-polarized waves
        % decouple. For cylinders this corresponds to decoupling of flexural and
        % longitudinal waves from torsional waves. 
        % Arguments: 
        % - n : circumferential wavenumber (only necessary for cylinders)
        % 
        % See also: Lamb, sh, fullyCoupled.
        lays = obj.lay;
        for l = lays % test each of the layers
            if ~l.decouplesLambvsSH(n)
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
        if ~isempty(obj.geom.gdofDBC)
            error('GEWTOOL:Waveguide:multipleCallsToFixGdof', ['Only one call to fixGdof is permitted.'... 
                ' Collect all global dofs you want to fix in a vector and pass them to fixGdof().']);
        end
        ops = fields(obj.op);
        for i=1:length(ops)
            opName = ops{i};
            obj.op.(opName)(gdof,:) = []; % remove row
            obj.op.(opName)(:,gdof) = []; % remove column
        end
        obj.geom.gdofDBC = gdof(:).'; 
    end

    function wc = cutoffFreq(obj, wmax, includeZeroFreq)
        if isempty(obj.op)
            error('GEWTOOL:cutoffFreq:chooseWaves', 'First choose the waves you want to compute: Lamb(), sh(), fullyCoupled().');
        end
        if size(obj.op.M,1) < 4
            error('GEWTOOL:cutoffFreq:chooseWaves', 'Matrices need to be at least of size 4x4.');
        end
        if nargin < 3
            includeZeroFreq = false;
        elseif ischar(includeZeroFreq)
            if strcmp(includeZeroFreq, 'includeZeroFreq')
                includeZeroFreq = true; 
            else
                includeZeroFreq = false;
            end
        elseif ~islogical(includeZeroFreq)
            error('GEWTOOL:cutoffFreq:wrongArg', 'Provide argument "includeZeroFreq" as a boolean/logical or char.');
        end
        wnCutoff = real(sqrt(eig(-obj.op.L0, obj.op.M, 'chol'))); % calculate cutoff frequencies
        wnCutoff = sort(wnCutoff);
        if nargin < 3 || ~includeZeroFreq % remove 0 Hz: relative to 4th cutoff as we have at most 3 cutoffs at 0 Hz.
            wnCutoff = wnCutoff(wnCutoff > 1e-4*wnCutoff(4));
        end
        wc = wnCutoff*obj.np.fh0/obj.np.h0; % in Hz.
        if nargin > 1 && ~isempty(wmax)
            wc = wc(wc <= wmax); % limit to maximum frequency if provided
        end
    end

    function nM = nModes(obj, wmax)
        nM = length(obj.cutoffFreq(wmax, 'includeZeroFreq'));
    end

    % declare functions implemented in external files:
	gew = assembleLayers(obj, udof, n)
end

end % class
