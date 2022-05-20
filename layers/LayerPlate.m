classdef LayerPlate < Layer
    
    methods
        function obj = LayerPlate(mat, ys, N)
            % LayerPlate: constructor
            obj = obj@Layer(mat, ys, N)
        end
        
        function [L0, L1, L2] = stiffnessOp(obj, udof, varargin)
            % stiffnessOp stiffness operator 
            cn = obj.mat.c/obj.mat.c(1,2,1,2); % normalized stiffness tensor
            % relevant material matrices: 
            cxx = squeeze(cn(1,udof,udof,1));
            cxy = squeeze(cn(1,udof,udof,2)); 
            cyx = squeeze(cn(2,udof,udof,1));
            cyy = squeeze(cn(2,udof,udof,2));
            % normalized element stiffness and flux:
            k2 = obj.PP; k1 = obj.PPd; % element stiffness
            g1 = -obj.PPd.'; g0 = obj.PdPd; % element flux
            % assemble:
            K2 = kron(cxx, k2); K1 = kron(cxy, k1); 
            G1 = kron(cyx, g1); G0 = kron(cyy, g0);
            % combine to polynomial of (ik):
            L0 = G0; L1 = K1 + G1; L2 = K2;
        end
        
        function [G0, G1] = tractionOp(obj, udof, varargin)
            % tractionOp traction operator 
            % TODO what is this good for?
            cn = obj.mat.c/obj.mat.c(1,2,1,2); % normalized stiffness tensor
            % relevant material matrices: 
            cyx = squeeze(cn(2,udof,udof,1));
            cyy = squeeze(cn(2,udof,udof,2));
            % assemble:
            g1 = -obj.PPd.'; g0 = obj.PdPd; % normalized element flux
            G1 = kron(cyx, g1); G0 = kron(cyy, g0); % assemble
        end
        
    end % methods
end % class
