classdef LayerPlate
    properties (Access = public)
        y   % collocation points on domain of unit length
        eta % collocation points on [-1, 1] TODO: is this being used?
        h   % thickness in m
        N   % number of collocation points
        D1  % diff matrix on unit domain 
        D2  % second order diff matrix
        mat % material 
    end
    
    methods
        function obj = LayerPlate(mat, h, N)
            % LayerPlate: constructor
            obj.mat = mat;
            obj.N = N;
            [obj.eta, D_dash] = chebdif(obj.N, 2);
            obj.D1 = 2*D_dash(:,:,1);
            obj.D2 = 4*D_dash(:,:,2);
            obj.h = h(2) - h(1);
            obj.y = obj.eta/2;
        end
        
        function [L0, L1, L2] = stiffnessOp(obj, udof, n)
            % stiffnessOp stiffness operator 
            cn = obj.mat.tensor/obj.mat.tensor(1,2,1,2); % normalized stiffness tensor
            warning('where is normalized stiffness tensor stored?')
            % relevant material matrices: 
            cxx = squeeze(cn(1,udof,udof,1));
            cxy = squeeze(cn(1,udof,udof,2)); 
            cyx = squeeze(cn(2,udof,udof,1));
            cyy = squeeze(cn(2,udof,udof,2));
            I = eye(size(cxx)); 
            % differentiation matrices on normalized domain:
            D1 = obj.D1; D2 = obj.D2; 
            Id = eye(size(D1)); % identity matrix for discretization
            % operators:
            L2 = kron(cxx, Id); L1 = kron(cxy + cyx, D1); L0 = kron(cyy, D2); 
        end
        
        function M = massOp(obj, udof)
            % massOp mass operator 
            rhon = 1*eye(length(udof)); % normalized mass matrix tensor 
            Id = eye(size(obj.D1)); % identity matrix for discretization
            M = kron(rhon, Id);
            warning('maybe put to super-class?')
        end
        
        function [B0, B1] = tractionOp(obj, udof, n)
            % tractionOp traction operator 
            cn = obj.mat.tensor/obj.mat.tensor(1,2,1,2); % normalized stiffness tensor
            % relevant material matrices: 
            cyx = squeeze(cn(2,udof,udof,1));
            cyy = squeeze(cn(2,udof,udof,2));
            % differentiation matrices on normalized domain:
            D1 = obj.D1; 
            Id = eye(size(obj.D1)); % identity matrix for discretization
            % operators:
            B1 = kron(cyx, Id([1, obj.N], :)); 
            B0 = kron(cyy, D1([1, obj.N], :));
        end
        
    end % methods
end % class