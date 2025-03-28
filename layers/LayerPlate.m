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
            K1 = kron(cxz, obj.PD); % stiffness
            G1 = kron(czx, -obj.PD.');  % boundary flux
            G0 = kron(czz, -obj.DD.'); % boundary flux
            % combine to polynomial of (ik):
            L0 = G0/hl;  L1 = K1 + G1;  L2 = K2*hl;
        end

        function decoupl = decouplesPolarization(obj,dof,~)
            % decouplesPolarization - Tests whether 'dof' decouples from the remaining 
            % degrees of freedom.
            %
            % Arguments: 
            % - dof : vector of degrees of freedom, e.g., [1 3] for Lamb waves.
            % - n   : circumferential wavenumber (only necessary for axial waves in cyl.)
            % 
            % See also: Waveguide.decouplesLambvsSH, Waveguide.decouplesPolarization
            rem = setdiff(1:3, dof); % remaining degrees of freedom 
            c1test = obj.mat.c(dof,dof,rem,dof); 
            c2test = obj.mat.c(dof,rem,dof,dof);
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
