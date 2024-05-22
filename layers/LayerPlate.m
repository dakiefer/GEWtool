classdef LayerPlate < Layer
% LayerPlate - Class to represent one layer of a multi-layered plate.
% There is usually no need to use this class explicitly (used internally by 
% Plate).
%
% See also Plate, Waveguide.
% 
% 2022 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France
    
    methods
        function obj = LayerPlate(mat, ys, N)
            % LayerPlate - Create a LayerPlate object.
            obj = obj@Layer(mat, ys, N);
        end
        
        function [L0, L1, L2] = stiffnessOp(obj, udof, hl, ~)
            % stiffnessOp - stiffness operators L0, L1, L2
            cn = obj.mat.c/obj.mat.c(1,2,1,2); % normalized stiffness tensor
            % relevant material matrices: 
            cxx = squeeze(cn(1,udof,udof,1));
            cxy = squeeze(cn(1,udof,udof,2)); 
            cyx = squeeze(cn(2,udof,udof,1)); % boundary flux
            cyy = squeeze(cn(2,udof,udof,2)); % boundary flux
            % assemble element stiffness terms:
            K2 = kron(cxx, obj.PP); 
            K1 = kron(cxy, obj.PPd); % stiffness
            G1 = kron(cyx, -obj.PPd.');  % boundary flux
            G0 = kron(cyy, -obj.PdPd.'); % boundary flux
            % combine to polynomial of (ik):
            L0 = G0/hl;  L1 = K1 + G1;  L2 = K2*hl;
        end
        
        function decoupl = decouplesLambvsSH(obj,~)
            % decouplesLambvsSH - Tests whether the Lamb- and SH-polarized waves decouple. Argument n is optional and does nothing.
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
