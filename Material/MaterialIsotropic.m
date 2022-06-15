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
        function mu = get.mu(obj)
            mu = obj.C(4,4);
        end
        function E = get.E(obj)
            [E, ~] = MaterialIsotropic.lame2Enu(obj.lambda, obj.mu);
        end
        function nu = get.nu(obj)
            [~, nu] = MaterialIsotropic.lame2Enu(obj.lambda, obj.mu);
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
    end
end