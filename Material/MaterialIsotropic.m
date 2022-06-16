classdef MaterialIsotropic < Material
    %MaterialIsotropic specialization for isotropic materials.

    properties (Dependent)
        lambda      % first Lamé parameter
        mu          % second Lamé parameter
        E           % Young's modulus
        nu          % Poisson's ratio
    end

    methods
        function obj = MaterialIsotropic(varargin)
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
        % Enu2lame Convert Young's modulus and Poisson ratio to Lamé parameters lambda
        % and mu.
            lbd = nu*E/((1+nu)*(1-2*nu));
            mu = E/(2*(1+nu));
        end
        function [E, nu] = lame2Enu(lbd, mu)
        % lame2Enu Convert Lamé parameters lambda and mu to Young's modulus E
        % and Poisson ratio nu.
            E = mu*(3*lbd + 2*mu)/(lbd + mu);
            nu = lbd/(2*(lbd + mu));
        end
        function c = lame2stiffnessTensor(lbd, mu)
            II = eye(3).*shiftdim(eye(3), -2); % 4th order "unit tensor"
            c = lbd*II + mu*(permute(II, [1 3 4 2]) + permute(II, [1 3 2 4])); % stiffness tensor
        end
    end
end