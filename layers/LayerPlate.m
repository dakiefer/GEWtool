classdef LayerPlate < Layer
% LayerPlate - Class to represent one layer of a multi-layered plate.
% There is usually no need to use this class explicitly (used internally by 
% Plate).
%
% See also Plate, Waveguide.
% 
% 2022-2024 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France
    
    methods
        function obj = LayerPlate(mat, zs, N)
            % LayerPlate - Create a LayerPlate object.
            obj = obj@Layer(mat, zs, N);
        end
        
        function [L0, L1, L2] = stiffnessOp(obj, udof, np, hl, ~)
            % stiffnessOp - stiffness operators L0, L1, L2
            cn = obj.mat.c/np.c0; % normalized stiffness tensor
            hl = hl/np.h0; % normalize
            % relevant material matrices: 
            cxx = squeeze(cn(1,udof,udof,1));
            cxz = squeeze(cn(1,udof,udof,3)); 
            czx = squeeze(cn(3,udof,udof,1)); % boundary flux
            czz = squeeze(cn(3,udof,udof,3)); % boundary flux
            % assemble element stiffness terms:
            K2 = kron(cxx, obj.PP); 
            K1 = kron(cxz, obj.PPd); % stiffness
            G1 = kron(czx, -obj.PPd.');  % boundary flux
            G0 = kron(czz, -obj.PdPd.'); % boundary flux
            % combine to polynomial of (ik):
            L0 = G0/hl;  L1 = K1 + G1;  L2 = K2*hl;
        end
        
        function decoupl = decouplesLambvsSH(obj,~)
            % decouplesLambvsSH - Tests whether the Lamb- and SH-polarized waves decouple. Argument n is optional and does nothing.
            lb = Waveguide.udofLamb;
            sh =  Waveguide.udofSH;
            c1test = obj.mat.c(lb,lb,sh,lb); 
            c2test = obj.mat.c(lb,sh,lb,lb);
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
