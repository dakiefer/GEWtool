classdef MaterialPiezoelectric < Material
% MaterialPiezoelectric - Represent piezoelectric material data. 
% Subclass of "Material". It represents material in the same way but
% additionally enables to access and manipulate the electrostatic and
% piezoelectric material parameters.
% 
% Example:
% mat = MaterialPiezoelectric('lithium_niobate')   % load from lithium_niobate.json (anywhere on path)
%
% See also: Material.
%
% 2024 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France

properties
    epsilon  % electric permittivity tensor
    eVoigt   % piezoelectric coupling parameters in Voigt notation
    eps0 = 8.8541878188e-12; % permittivity of vacuum
end

properties (Dependent)
    e        % piezoelectric coupling parameters (3rd order tensor)
end

methods
    function obj = MaterialPiezoelectric(varargin)
        % MaterialPiezoelectric - Create a piezoelectric material object.
        % 
        % Usage:
        % mat = MaterialPiezoelectric('lithium_niobate');   % load from file "lithium_niobate.json" (somewhere on path)
        % mat = MaterialPiezoelectric(param); % from structure "param" with fields "name", "C", "rho", "eVoigt", and "epsl"
        %
        if isa(varargin{1}, 'Material') || isa(varargin{1}, 'struct') % convert from subclasses to superclass
            data = varargin{1};
        elseif ischar(varargin{1}) || isstring(varargin{1})  % load by name
            matname = char(varargin{1});
            filename = [matname '.json'];
            data = jsondecode(fileread(filename));
        else
            error('GEWTOOL:Material:unknownArgument', 'Provide filename to load or structure describing material.');
        end
        obj = obj@Material(data);
        obj.eVoigt = data.e;
        obj.epsilon = data.epsr*obj.eps0;
    end

    function matStruct = struct(obj) 
        % struct - convert Material object to a struct.
        matStruct.name = obj.name;
        matStruct.rho = obj.rho;
        matStruct.C = obj.C;
        matStruct.e = obj.eVoigt; 
        matStruct.epsr = obj.epsilon/obj.eps0; 
        matStruct.symmetry = obj.symmetry;
    end

    function e = get.e(obj)
        e = voigt2tensor(obj.eVoigt);
    end
    function obj = set.e(obj, e)
        validateattributes(e,'numeric',{'size',[3 3 3]});
        obj.eVoigt = tensor2voigt(e);
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
        obj = permute@Material(obj, perm);
        obj.epsilon = obj.epsilon(perm,perm);
        obj.e = obj.e(perm,perm,perm);
    end

    function obj = transformBasis(obj, Q)
        % TRANSFORMBASIS - Change of basis "Q" applied to all tensor quantities.
        % Symmetries are conserved and numerical errors are removed by rounding
        % to 12 digits with respect to the Frobenius norm.
        %
        % Arguments:
        % - obj: Material object.
        % - Q:   Orthogonal matrix discribing the transformation, e.g., a rotation.
        c0 = norm(obj.c(:)); % normalization constant
        crot = transformBasis(obj.c/c0, Q);
        crot = round( crot, floor(-log10(eps))-3 )*c0; % 12 digits accuracy + upscaling
        e0 = norm(obj.e(:));
        erot = transformBasis(obj.e/e0, Q);
        erot = round( erot, floor(-log10(eps))-3 )*e0; % 12 digits accuracy + upscaling
        E0 = norm(obj.epsilon(:));
        Erot = transformBasis(obj.epsilon/E0, Q); 
        Erot = round (Erot, floor(-log10(eps))-3 )*E0; 
        if obj.isSymmetric
            crot = (crot + permute(crot, [3 4 1 2]))/2; % avoid assymetry due to rotation
        end
        if obj.isSymmetric([2 1 3 4])
            crot = (crot + permute(crot, [2 1 3 4]))/2; % avoid assymetry due to rotation
        end
        if obj.isSymmetric([1 2 4 3])
            crot = (crot + permute(crot, [1 2 4 3]))/2; % avoid assymetry due to rotation
        end
        if all( abs( obj.e - permute(obj.e,[1 3 2]) ) <= 1e4*eps*e0, 'all')
            erot = (erot + permute(erot, [1 3 2]))/2; % avoid assymetry due to rotation
        end
        if all( abs( obj.epsilon - obj.epsilon.' ) <= 1e4*eps*E0, 'all')
            Erot = (Erot + Erot.')/2; % avoid assymetry due to rotation
        end
        obj.c = crot;
        obj.e = erot; 
        obj.epsilon = Erot;
    end
    
    function sym = isInvariantOnReflection(obj, ax, tol)
        % ISINVARIANTONREFLECTION - test if invariant to reflection along axis 'ax'.
        % 
        % Arguments:
        % - obj:   Material object.
        % - ax:    axis as a scalar integer 1 = x, 2 = y, or 3 = z. (Default: 3)
        % - tol:   tolerance for finite precision test. (Default: 1e4*eps)
        if nargin < 3
            tol = 1e4*eps;
        end
        if nargin < 2
            ax = 3; % default is z-axis
        end
        matR = obj.reflect(ax);
        matOrig = obj.transformBasis(eye(3)); % same rounding as in matR
        csym = all((matOrig.c - matR.c)/norm(matOrig.c(:)) <= tol, 'all'); % all zero -> true
        Esym = all((matOrig.eps - matR.eps)/norm(matOrig.eps(:)) <= tol, 'all'); % all zero -> true
        esym = all((matOrig.e - matR.e)/norm(matOrig.e(:)) <= tol, 'all'); % all zero -> true
        sym = csym & Esym & esym;
    end

    function dis = isDissipative(obj)
        dis = true;
        cTest = isreal(obj.c) && obj.isSymmetric; 
        ETest = isreal(obj.epsilon) && issymmetric(obj.epsilon); 
        eTest = isreal(obj.e) && all( abs( obj.e - permute(obj.e,[1 3 2]) ) <= 1e4*eps*norm(obj.e(:)), 'all');
        if cTest && ETest && eTest
            dis = false; 
        end
    end

    function cons = isConservative(obj)
        cons = ~isDissipative(obj);
    end

    function save(obj,path)
        save@Material(obj,path);
    end

    % maye the following should also be specialized for piezoelectric media: 
    % [s, e0] = slownessCurve(obj, alpha, ek)
    % [cs, eu] = wavespeeds(obj, ek)
    % [ce] = energyVel(obj, ek)
    % plotSlownessCurve(obj, erot, varargin)
    % plotRayCurve(obj, erot, varargin)
    

    %% overload operators: 
    function ret = eq(a, b)
        % eq - Test if the material parameters of "a" and "b" are the same up to numerics.
        % Usage: 
        % isEq = eq(a, b);
        % isEq = a == b;
        ret = ~any(a.c ~= b.c, 'all') && ...
               a.rho == b.rho &&...
              ~any(a.e ~= b.e, 'all') && ...
              ~any(a.epsilon ~= b.epsilon, 'all');
    end


end % methods

end % class