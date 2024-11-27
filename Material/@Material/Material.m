classdef (InferiorClasses = {?MaterialIsotropic}) Material
% Material - Represent mechanical material data (generally anisotropic).
% Stores and manipulates elasticity moduli, density and derived quantities. The
% full 4th order stiffness tensor c is stored as a [3x3x3x3] array (helpful for
% nonlinear materials). It interfaces to a simple material database that allows
% to load data by name.
% 
% Conversion functions between sets of material parameters are provided by the
% subclass MaterialIsotropic.
% 
% Examples:
% mat = Material('steel')   % load from steel.json (anywhere on path)
% similarMat = Material('anyName', mat.C, 1.1*mat.rho);   % based on steel
% 
% See also: Material.Material, MaterialIsotropic.
%
% 2022-2024 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France

properties
    name        % name for the material: string
    symmetry    % symmetry class: e.g., "isotropic": string
    c double    % 4th order stiffness tensor
    rho double  % mass density [1x1]
end

properties (Dependent)
    C double % Voigt notated stiffness [6x6]
    cl       % (quasi)-longitudinal wave speed in x-direction
    ct       % first (quasi)-transverse wave speed in x-direction
    ct2      % second (quasi)-transverse waves speed in x-direction
    AU       % universal anisotropy index (doi 10.1103/PhysRevLett.101.055504)
end

methods
    function obj = Material(varargin)
        % Material - Create a Material object. 
        %
        % Usage:
        % mat = Material('steel');   % load from file "steel.json" (somewhere on path)
        % mat = Material('name', C, rho);  % from 6x6 Voigt stiffness and density
        % mat = Material('name', c, rho);  % from 3x3x3x3 stiffness and density
        % mat = Material('name', lbd, mu, rho);  % from lamé parameters and density
        % mat = Material(param); % from structure param with fields "name", "C" or "c", "rho"
        % mat = Material(mat); % from Material object "mat" (conversion from subclasses)
        %
        if nargin == 1
            if isa(varargin{1}, 'Material') || isa(varargin{1}, 'struct') % convert from subclasses to superclass
                data = varargin{1};
            elseif ischar(varargin{1}) || isstring(varargin{1})  % load by name
                matname = char(varargin{1});
                filename = [matname '.json'];
                data = jsondecode(fileread(filename));
            else
                error('GEWTOOL:Material:unknownArgument', 'Provide filename to load or structure describing material.');
            end
            % add missing fields:
            if isstruct(data) && ~isfield(data, 'symmetry') % add field default value
                if isfield(data, 'lambda') && isfield(data, 'mu')
                    data.symmetry = 'isotropic';
                else
                    data.symmetry = 'unknown'; % might still be isotropic... 
                end
            end
            if isstruct(data) && ~isfield(data, 'name') % add field default value
                data.name = 'unnamed';
            end
        elseif nargin == 3  % name, C or c, rho 
            validateattributes(varargin{1},{'char'},{'vector'},1); % name
            data.name = varargin{1};
            if ~isnumeric(varargin{2}) % stiffness
                error('Expected stiffness (argument 2) to be numeric.');
            elseif isequal(size(varargin{2}), [6,6]) % Voigt noteted tensor
                data.C = varargin{2};
            elseif isequal(size(varargin{2}), [3,3,3,3]) % 4th order tensor
                data.c = varargin{2};
            else
                error('Expected argument 2 to be of size [6x6] or [3x3x3x3]');
            end
            validateattributes(varargin{3},{'numeric'},{'size',[1 1]},3); % density
            data.rho = varargin{3};
            data.symmetry = 'anisotropic';
        elseif nargin == 4 && (ischar(varargin{1}) || isstring(varargin{1})) % name, lambda, mu and rho 
            validateattributes(varargin{2},{'numeric'},{'size',[1 1]});
            validateattributes(varargin{3},{'numeric'},{'size',[1 1]});
            validateattributes(varargin{4},{'numeric'},{'size',[1 1]});
            data.name = varargin{1}; 
            data.lambda = varargin{2};
            data.mu = varargin{3};
            data.rho = varargin{4};
            data.symmetry = 'isotropic';
        end
        if strcmp(data.symmetry, 'isotropic') 
            data.c = MaterialIsotropic.lame2stiffnessTensor(data.lambda, data.mu);
        elseif isfield(data, 'C') 
            data.c = voigt2tensor(data.C);
        end
        obj.c = data.c;
        obj.rho = data.rho;
        obj.name = data.name;
        obj.symmetry = data.symmetry;
        if any(diag(obj.C) <= 0)
            warning('GEWTOOL:Material:stability', 'The material has nonpositive diagonal entries.');
        elseif any(eig(obj.C) <= 0)
            warning('GEWTOOL:Material:stability', 'The material exhibits nonpositive eigenvalue(s).');
        end
    end

    function matStruct = struct(obj) 
        % struct - convert Material object to a struct.
        matStruct.name = obj.name;
        matStruct.rho = obj.rho;
        if strcmp(obj.symmetry, 'isotropic') 
            matStruct.lambda = obj.lambda;
            matStruct.mu = obj.mu;
        else
            matStruct.C = obj.C;
        end
        matStruct.symmetry = obj.symmetry;
    end

    function C = get.C(obj)
        C = tensor2voigt(obj.c);
    end
    function obj = set.C(obj, C)
        obj.c = voigt2tensor(C);
    end
    function cl = get.cl(obj)
        cs = wavespeeds(obj);
        cl = cs(1);
    end
    function ct = get.ct(obj)
        cs = wavespeeds(obj);
        ct = cs(2);
    end
    function ct2 = get.ct2(obj)
        cs = wavespeeds(obj);
        ct2 = cs(3);
    end
    function AU = get.AU(obj) 
        if obj.symmetry == "isotropic" || ...
                (obj.symmetry == "cubic" && length(unique(obj.C(obj.C ~= 0))) == 3) % ensures principal axis
            A = 2*obj.C(4,4)/(obj.C(1,1)-obj.C(1,2)); % Zener anisotropy index
            AU = 6/5*(sqrt(A) - 1/sqrt(A))^2;
        else
            AU = AUanisotropyIndex(obj,1200);
        end
        AU = round(AU,3); 
    end

    function obj = permute(obj, perm)
        % permute - Permute the material coordinate system. 
        % The default transformation is ex-ey-ez -> ez-ex-ey, i.e. 
        % 1->3, 2->1, 3->2, given by perm = [3,1,2].
        %
        % Argument:
        % - perm:   [3x1]-vector indicating the permutation 1:3 -> perm.
        if nargin < 2
            perm = [3,1,2]; % 3->1, 1->2, 2->3
        end
        obj.c = obj.c(perm,perm,perm,perm);
    end

    function sym = isSymmetric(obj, perm, tol)
        % ISSYMMETRICONPERM - Test symmetry upon tensor index permutation.
        % 
        % Usage: 
        % sym = obj.isSymmetric();       % default perm = [3,4,1,2]
        % sym = obj.isSymmetric(perm);
        % sym = obj.isSymmetric(perm, tol); % also provide rel. tolerance for testing
        % 
        % Argument: 
        % - perm: [4x1]-vector of stiffness tensor index permuation 1:4 -> perm
        %         default: perm = [3,4,1,2] (i.e., cijkl -> cklij).
        % - tol:  (scalar, default: 0) relative tolerance for testing symmetry
        if nargin < 3
            tol = 0;
        end
        if nargin < 2 || isempty(perm)
            perm = [3,4,1,2];
        end
        cperm = permute(obj.c, perm);
        sym = all( abs(cperm - obj.c) <= tol*norm(obj.c(:)), 'all');
    end

    function obj = rotateEuler(obj, varargin)
        % ROTATEEULER - Rotate the stiffness tensor by the given angle-axis sequence.
        %
        % Arguments:
        % - obj:     Material object.
        % - angle:   (scalar numeric) angle in radian to rotate
        % - axis:    (one of 'x', 'y', 'z') axis to rotate about
        % You can specify as many sets of angle-axis pairs as you desire. The rotation
        % will be performed in the provided order in the fixed initial coordinate system
        % (extrinsic).
        % 
        % Alternative Arguments:
        % - A:       (3x3x...x3 numeric) nth-order tensor (all dimensions 3)
        % - angleX:  (scalar numeric) angle in rad to turn around axis x
        % - angleY:  (scalar numeric) angle in rad to turn around axis y
        % - angleZ:  (scalar numeric) angle in rad to turn around axis z
        % This is an extrinsic, passive rotation around x-y-z (in that order). Thereby, 
        % x-y-z is the original fixed coordinate system. Note that these are improper
        % Euler angles (also called Cardan angles or Tait–Bryan angles).
        % For more details see: https://en.wikipedia.org/wiki/Euler_angles
        % 
        % usage: 
        % R = rotateEuler(A, angle, axis); 
        % R = rotateEuler(A, angle1, axis1, angle2, axis2); % rotate around
        %          "axis1", then around "axis2".
        % R = rotateEuler(A, angleX, angleY, angleZ); % rotate x, then y, then z.
        % 
        % Literature: D. Royer and T. Valier-Brasier, Ondes élastiques dans les solides 
        % (Elastic waves in solids), vol. 1, 2 vols. London: ISTE éditions, 2021.
        %
        % 2022-2024 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France
        R = eulerAnglesToRotationMatrix(varargin{:});
        obj = transformBasis(obj,R); % 
    end

    function obj = reflect(obj, ax)
        % REFLECT - Reflection of material along axis 'ax'.
        %
        % Arguments:
        % - obj:   Material object
        % - ax:    axis as a scalar integer 1 = x, 2 = y, or 3 = z.
        Qx = [-1,  0,   0;
              0,   1,   0;
              0,   0,   1]; % reflexion ex -> -ex
        Q = circshift(Qx, [ax-1, ax-1]);
        obj = obj.transformBasis(Q); % keeps symmetry
    end

    function obj = transformBasis(obj, Q)
        % TRANSFORMBASIS - Change of basis "Q" applied to stiffness tensor.
        % Symmetries are conserved and numerical errors are removed by rounding
        % to 12 digits with respect to the Frobenius norm.
        %
        % Arguments:
        % - obj: Material object.
        % - Q:   Orthogonal matrix discribing the transformation, e.g., a rotation.
        cn = norm(obj.c(:)); % normalization constant
        crot = transformBasis(obj.c/cn, Q);
        crot = round( crot, floor(-log10(eps))-3 )*cn; % 12 digits accuracy + upscaling
        if obj.isSymmetric
            crot = (crot + permute(crot, [3 4 1 2]))/2; % avoid assymetry due to rotation
        end
        if obj.isSymmetric([2 1 3 4])
            crot = (crot + permute(crot, [2 1 3 4]))/2; % avoid assymetry due to rotation
        end
        if obj.isSymmetric([1 2 4 3])
            crot = (crot + permute(crot, [1 2 4 3]))/2; % avoid assymetry due to rotation
        end
        obj.c = crot;
    end
    
    function sym = isInvariantOnReflection(obj, ax, tol)
        % ISINVARIANTONREFLECTION - test if invariant to reflection along axis 'ax'.
        % 
        % Arguments:
        % - obj:   Material object.
        % - ax:    axis as a scalar integer 1 = x, 2 = y, or 3 = z. (Default: 2)
        % - tol:   tolerance for finite precision test. (Default: 1e4*eps)
        if nargin < 3
            tol = 1e4*eps;
        end
        if nargin < 2
            ax = 3; % default is z-axis
        end
        matR = obj.reflect(ax);
        matOrig = obj.transformBasis(eye(3)); % same rounding as in matR
        sym = all((matOrig.c - matR.c)/matOrig.c(1,1,1,1) <= tol, 'all'); % all zero -> true
    end

    function decoupl = decouplesSA(obj)
        % DECOUPLESSA - test if invariant to reflection along y-axis.
        % Basically an alias to isInvariantOnReflection().
        decoupl = obj.isInvariantOnReflection(3); % full polarization
    end

    function dis = isDissipative(obj)
        dis = true;
        if isreal(obj.c) && obj.isSymmetric
            dis = false;
        end
    end

    function cons = isConservative(obj)
        cons = false;
        if isreal(obj.c) && obj.isSymmetric
            cons = true;
        end
    end

    function cmin = minWavespeed(obj, en)
        % CMIN - Approximate minimum plane wave speed in the en-plane. 
        if nargin < 2
            en = [0; 0; 1];
        end
        alpha = linspace(0, pi/2, 100); % TODO this might be assuming some symmetries in the material 
        s = obj.slownessCurve(alpha, en);
        cmin = 1/max(s(:));
    end

    function cmax = maxWavespeed(obj, en)
        % CMIN - Approximate maximum plane wave speed in the en-plane. 
        if nargin < 2
            en = [0; 0; 1];
        end
        alpha = linspace(0, pi/2, 100); % TODO this might be assuming some symmetries in the material 
        s = obj.slownessCurve(alpha, en);
        cmax = 1/min(s(:));
    end

    [s, e0] = slownessCurve(obj, alpha, ek)
    [cs, eu] = wavespeeds(obj, ek)
    [ce] = energyVel(obj, ek)
    plotSlownessCurve(obj, erot, varargin)
    plotRayCurve(obj, erot, varargin)
    save(obj,path)

    %% overload operators: 
    function ret = eq(a, b)
        % eq - Test if stiffness c and density rho are the same for materials a and b.
        % Usage: 
        % isEq = eq(a, b);
        % isEq = a == b;
        ret = ~any(a.c ~= b.c, 'all') && a.rho == b.rho;
    end
    function ret = ne(a, b)
        ret = ~eq(a, b);
    end
    function plot(varargin)
        % plot - plot the slowness curve around the axis ex.
        % Alias to plotSlownessCurve(varargin{:});
        plotSlownessCurve(varargin{:});
    end
    
end % methods
end % class
