classdef LayerPlate < Layer
% LayerPlate - Class to represent one layer of a multi-layered plate.
% There is usually no need to use this class explicitly (used internally by 
% Plate).
%
% See also Plate, Waveguide.
% 
% 2022 - Daniel A. Kiefer
% Institut Langevin, Paris, France
% 
    
    methods
        function obj = LayerPlate(mat, ys, N)
            % LayerPlate - Create a LayerPlate object.
            obj = obj@Layer(mat, ys, N);
        end
        
        function [L0, L1, L2] = stiffnessOp(obj, udof, varargin)
            % stiffnessOp - stiffness operators L0, L1, L2
            cn = obj.mat.c/obj.mat.c(1,2,1,2); % normalized stiffness tensor
            % relevant material matrices: 
            cxx = squeeze(cn(1,udof,udof,1));
            cxy = squeeze(cn(1,udof,udof,2)); 
            % assemble:
            K2 = kron(cxx, obj.PP/obj.h); K1 = kron(cxy, obj.PPd)/obj.h; % stiffness
            [G0, G1] = obj.tractionOp(udof, varargin{:}); % flux
            % combine to polynomial of (ik):
            L0 = G0; L1 = K1 + G1; L2 = K2;
        end
        
        function [G0, G1] = tractionOp(obj, udof, varargin)
            % tractionOp - traction operator (flux, used internally)
            cn = obj.mat.c/obj.mat.c(1,2,1,2); % normalized stiffness tensor
            % relevant material matrices: 
            cyx = squeeze(cn(2,udof,udof,1));
            cyy = squeeze(cn(2,udof,udof,2));
            % normalized element flux
            G1 = kron(cyx, -obj.PPd.')/obj.h; G0 = kron(cyy, -obj.PdPd.')/obj.h; % assemble
        end

        function M = massOp(obj, udof)
            % massOp - mass operator M
            rhon = eye(length(udof)); % normalized mass matrix (for each dof in u) 
            me = obj.PP; % element mass
            M = kron(rhon,me)/obj.h; % assemble
        end
        
        function decoupl = decouplesLambvsSH(obj)
            % decouplesLambvsSH - Tests whether the Lamb- and SH-polarized waves decouple.
            xy = [1 2]; % Lamb polarization
            z =  3; % SH polarization 
            c1test = obj.mat.c(xy,xy,z,xy); 
            c2test = obj.mat.c(xy,z,xy,xy);
            decoupl = all(c1test(:) == 0) & all(c2test(:) == 0);
        end

        %% overload operators: 
        function ret = eq(a, b)
            ret = eq@Layer(a, b);
        end
        function ret = ne(a, b)
            ret = ne@Layer(a, b);
        end
        
    end % methods
end % class
