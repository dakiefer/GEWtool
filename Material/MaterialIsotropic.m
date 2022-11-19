classdef MaterialIsotropic < Material
% MaterialIsotropic - Represent isotropic mechanical material data.
% Subclass of "Material". It represents material in the same way but enables
% to access and manipulate the common isotropic elasticity parameters.
% Also provides static methods for conversion between different set of material
% parameters.
% 
% Example:
% mat = MaterialIsotropic('steel')   % load from steel.json (anywhere on path)
% heavierMat = MaterialIsotropic('anyName', mat.lambda, mat.mu, 1.1*mat.rho);
%
% See also: MaterialIsotropic.MaterialIsotropic, Material.
%
% 2022 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France

    properties (Dependent)
        lambda      % first Lamé parameter
        mu          % second Lamé parameter
        E           % Young's modulus
        nu          % Poisson's ratio
    end

    methods
        function obj = MaterialIsotropic(varargin)
            % MATERIALISOTROPIC - Create an isotropic material object.
            % 
            % Usage:
            % mat = MaterialIsotropic('steel');   % load from file "steel.json" (somewhere on path)
            % mat = MaterialIsotropic('name', lbd, mu, rho);  % from lamé parameters and density
            % mat = MaterialIsotropic(param); % from structure param with fields "name", "lambda", "mu" and "rho"
            % mat = MaterialIsotropic(mat); % from Material object "mat" (conversion from subclasses)
            %
            
            obj = obj@Material(varargin{:});
            if ~strcmp(obj.symmetry, 'isotropic') 
                error('Given material data is anisotrpic.');
            end
        end
        function lambda = get.lambda(obj)
            lambda = obj.C(1,2);
        end
        function obj = set.lambda(obj, lbd)
            obj.c = MaterialIsotropic.lame2stiffnessTensor(lbd, obj.mu);
        end
        function mu = get.mu(obj)
            mu = obj.C(4,4);
        end
        function obj = set.mu(obj, mu) 
            obj.c = MaterialIsotropic.lame2stiffnessTensor(obj.lambda, mu);
        end
        function E = get.E(obj)
            [E, ~] = MaterialIsotropic.lame2Enu(obj.lambda, obj.mu);
        end
        function obj = set.E(obj, E)
            [lbd, mu] = MaterialIsotropic.Enu2lame(E, obj.nu);
            obj.c = MaterialIsotropic.lame2stiffnessTensor(lbd, mu);
        end
        function nu = get.nu(obj)
            [~, nu] = MaterialIsotropic.lame2Enu(obj.lambda, obj.mu);
        end
        function obj = set.nu(obj, nu)
            [lbd, mu] = MaterialIsotropic.Enu2lame(obj.E, nu);
            obj.c = MaterialIsotropic.lame2stiffnessTensor(lbd, mu);
        end
    end

    methods (Static)
        function [lbd, mu] = Enu2lame(E, nu)
        % Enu2lame - Convert Young's modulus E and Poisson ratio nu to Lamé parameters lambda
        % and mu.
        % Input: scalars E and nu
        % Output: scalars lbd, mu
            lbd = nu*E/((1+nu)*(1-2*nu));
            mu = E/(2*(1+nu));
        end
        function [E, nu] = lame2Enu(lbd, mu)
        % lame2Enu - Convert Lamé parameters lambda and mu to Young's modulus E
        % and Poisson ratio nu.
        % Input: scalars lbd, mu
        % Output: scalars E, nu
            E = mu*(3*lbd + 2*mu)/(lbd + mu);
            nu = lbd/(2*(lbd + mu));
        end
        function c = lame2stiffnessTensor(lbd, mu)
            % lame2stiffnessTensor - Convert Lamé parameters lambda and mu to c 
            % (4th order stiffness tensor).
            % Input: scalars lbd, mu
            % Output: c [3x3x3x3]
            II = eye(3).*shiftdim(eye(3), -2); % 4th order "unit tensor"
            c = lbd*II + mu*(permute(II, [1 3 4 2]) + permute(II, [1 3 2 4])); % stiffness tensor
        end
    end
end