classdef LayerCylindrical 
    properties (Access = public)
        r   % collocation points on domain of unit length
        eta % collocation points on [-1, 1] TODO: is this being used?
        h   % thickness in m
        N   % number of collocation points
        D1  % diff matrix on unit domain 
        D2  % second order 
        mat % material 
        % n = 0 % circumferential order: exp(1i*n*phi)
    end

    methods
        function [obj] = LayerCylindrical(mat, rs, N)
            %% LayerSolid: constructor method
            obj.mat = mat;
            
            obj.N = N;
            [obj.eta, D_dash] = chebdif(obj.N, 2);
            obj.D1 = 2*D_dash(:,:,1);
            obj.D2 = 4*D_dash(:,:,2);
            
            obj.h = rs(2) - rs(1);
            obj.r = (obj.eta + 1)/2 + rs(1)/obj.h;
        end

        function [L0, L1, L2] = stiffnessOp(obj, dof, n)
            % stiffnessOp stiffness operator describing the left side of the 
            % elastodynamic wave equation:
            % L u = (w h0)^2 M u, 
            % where L = (i k_x h0)^2 L_2 + (i k_x h) L_1 + L_0$.

            cn = obj.mat.tensor/obj.mat.tensor(1,2,1,2); % normalized stiffness tensor
            % relevant material matrices: 
            crr = squeeze(cn(1,dof,dof,1));
            cpp = squeeze(cn(2,dof,dof,2));
            czz = squeeze(cn(3,dof,dof,3));
            Czp = squeeze(cn(3,dof,dof,2)) + squeeze(cn(2,dof,dof,3));
            Crz = squeeze(cn(3,dof,dof,1)) + squeeze(cn(1,dof,dof,3));
            Crp = squeeze(cn(1,dof,dof,2)) + squeeze(cn(2,dof,dof,1));
            crp = squeeze(cn(1,dof,dof,2));
            crz = squeeze(cn(1,dof,dof,3));
            % differetiation in curvilinear coordinate system:
            A = [0, -1, 0; 1, 0, 0; 0, 0, 0]; A = squeeze(A(dof, dof));
            B = [1,  0, 0; 0, 1, 0; 0, 0, 0]; B = squeeze(B(dof, dof));
            % differentiation matrices on normalized domain 
            Dr1 = obj.D1; Dr2 = obj.D2; r = obj.r;
            Id = eye(size(Dr1));                % identity matrix for discretization
            %% problem setup: (i*kh)^2*L2 + (i*kh)*L1(i*n) + L0(i*n) + w^2*M = 0
            L2 = kron(czz, Id); 
            L1 = kron(Czp*A, diag(1./r)) + kron(Crz, Dr1) + 1i*n*kron(Czp, diag(1./r));
            L0 = kron(crr, Dr2) + kron(Crp*A, diag(1./r)*Dr1 - diag(1./r.^2)) ...
                + kron(cpp, diag(1./r)*Dr1) - kron(cpp*B, diag(1./r.^2)) ...
                + 1i*n*(kron(Crp, diag(1./r)*Dr1 - diag(1./r.^2)) + kron(2*cpp*A, diag(1./r.^2))) ...
                + (1i*n).^2*(kron(cpp, diag(1./r.^2)));
        end

        function [M] = massOp(obj, dof)
            % massOp mass operator describing the right side of the elastodynamic
            % wave equation:
            % L u = -(w h0)^2 M u, 
            % where M = rho*d_ik. 
            
            rhon = 1*eye(length(dof)); % normalized mass matrix tensor 
            Id = eye(size(obj.D1));                % identity matrix for discretization
            M = kron(rhon, Id); 
        end

        function [B0, B1] = tractionOp(obj, dof, n)
            % tractionOp traction operator to impose traction boundary conditions.
            % The traction operator B is given by
            % B u = T ey,
            % with B = (i k_x h0)^2 B2 + (i k_x h0) B1 + B0,
            % where T is the stress tensor and ey is the unit vector normal to 
            % the plate.

            cn = obj.mat.tensor/obj.mat.tensor(1,2,1,2); % normalized stiffness tensor
            crr = squeeze(cn(1,dof,dof,1));
            crp = squeeze(cn(1,dof,dof,2));
            crz = squeeze(cn(1,dof,dof,3));
            A = [0, -1, 0; 1, 0, 0; 0, 0, 0]; A = squeeze(A(dof, dof)); % differetiation in curvilinear coordinate system
            % differentiation matrices on normalized domain 
            Dr1 = obj.D1; rInv = diag(1./obj.r); 
            Id = eye(obj.N);                % identity matrix for discretization
            % Assemble the mass operator M:
            B1 = kron(crz, Id([1, obj.N], :)); % BC going into L1
            B0 = kron(crr, Dr1([1, obj.N], :)) + kron(crp*A, rInv([1, obj.N], :)) + 1i*n*kron(crp, rInv([1, obj.N], :)); % BC going into L0
        end

        function [U0] = displacementOp(obj, dof)
            I = eye(length(dof));
            Id = eye(obj.N);  
            U0 = kron(I, Id([1, obj.N], :));
        end
        
    end
end
