classdef LayerPlate < Layer
    
    methods
        function obj = LayerPlate(mat, ys, N)
            % LayerPlate: constructor
            obj = obj@Layer(mat, ys, N)
        end
        
        function [L0, L1, L2] = stiffnessOp(obj, udof, varargin)
            % stiffnessOp stiffness operator 
            cn = obj.mat.c/obj.mat.c(1,2,1,2); % normalized stiffness tensor
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
        
        function [B0, B1] = tractionOp(obj, udof, varargin)
            % tractionOp traction operator 
            cn = obj.mat.c/obj.mat.c(1,2,1,2); % normalized stiffness tensor
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
