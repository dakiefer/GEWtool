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
            E = obj.mu*(3*obj.lambda + 2*obj.mu)/(obj.lambda + obj.mu);
        end
        function nu = get.nu(obj)
            nu = obj.lambda/(2*(obj.lambda + obj.mu));
        end
    end
end