classdef LayerPlatePiezo < LayerPlate
% LayerPlatePiezo - Class to represent one layer of a multi-layered piezoelectric plate.
% There is usually no need to use this class explicitly (used internally by 
% Plate). If the material "mat" is piezoelectric, the Plate class will
% automatically choose to construct a "LayerPlatePiezo". 
%
% See also Plate, Waveguide.
% 
% 2024-2025 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France
    
    methods
        function obj = LayerPlatePiezo(mat, zs, N)
            % LayerPlate - Create a LayerPlate object.
            if ~isa(mat, 'MaterialPiezoelectric')
                error('GEWTOOL:LayerPlatePiezo', 'The material should be of class MaterialPiezoelectric.');
            end
            obj = obj@LayerPlate(mat, zs, N);
        end

        function [M] = massOp(obj, udof, np, hl, ~)
            % assemble mass for all equations: 
            rhon = obj.mat.rho/np.rho0;
            MM  = blkdiag(rhon*eye(length(udof)), 0);  % mechanics only, singular
            M = kron(MM,obj.PP)*hl/np.h0; % assemble
        end
        
        function [L0, L1, L2] = stiffnessOp(obj, udof, np, hl, ~)
            % stiffnessOp - stiffness operators L0, L1, L2
            cn = obj.mat.c/np.c0; % normalized stiffness tensor
            en = obj.mat.e/np.e0;
            pn = obj.mat.epsilon/np.eps0; % permittivity: np.eps0 is average of epsilon tensor components
            hl = hl/np.h0; % normalize
            % relevant material matrices: 
            cxx = squeeze(cn(1,udof,udof,1));
            cxz = squeeze(cn(1,udof,udof,3));
            czx = squeeze(cn(3,udof,udof,1));
            czz = squeeze(cn(3,udof,udof,3));
            exx = squeeze(en(1,1,udof));
            exz = squeeze(en(1,3,udof));
            ezx = squeeze(en(3,1,udof));
            ezz = squeeze(en(3,3,udof));
            pxx = pn(1,1);
            pxz = pn(1,3);
            pzx = pn(3,1);
            pzz = pn(3,3);
            % assemble equations:
            % equation 1: balance of linear momentum, i.e., div(c:grad(u) + grad(Phi).e) + w^2 rho u = 0
            % equation 2: charge-free material, i.e.,       div(e:grad(u) - eps.grad(Phi)) = 0
            A  = [cxx,    exx;      % ~ (ik)^2
                  exx.'  -pxx];
            B  = [czx,    exz;      % ~ (ik)^1
                  ezx.'  -pzx];
            Bt = [cxz,    ezx;      % ~ (ik)^1
                  exz.'  -pxz];
            C  = [czz,    ezz;      % ~ (ik)^0
                  ezz.'  -pzz];
            
            % assemble element stiffness terms:
            K2 = kron(A,  obj.PP); 
            K1 = kron(Bt, obj.PD);
            G1 = kron(B, -obj.PD.'); % boundary flux
            G0 = kron(C, -obj.DD.'); % boundary flux
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

            if ~decouplesPolarization@LayerPlate(obj,dof)
                decoupl = false; return;
            end
            rem = setdiff(1:3, dof); % remaining degrees of freedom 
            etest = obj.mat.e(dof,rem,dof);
            decoupl = all(etest(:) == 0);
        end
        
    end % methods

    methods (Static) 
        function N = Nunknowns(udof)
            N = length(udof) + 1; 
        end
    end
end % class
