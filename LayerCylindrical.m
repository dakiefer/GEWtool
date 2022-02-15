classdef LayerCylindrical 
    properties (Access = public)
        r   % collocation points on domain of unit length
        eta % collocation points on [-1, 1] TODO: is this being used?
        h   % thickness in m
        N   % number of collocation points
        D1  % diff matrix on unit domain 
        D2  % second order diff matrix
        mat % material 
    end

    methods
        function obj = LayerCylindrical(mat, rs, N)
            % LayerCylindrical: constructor
            obj.mat = mat;
            obj.N = N;
            [obj.eta, D_dash] = chebdif(obj.N, 2);
            obj.D1 = 2*D_dash(:,:,1);
            obj.D2 = 4*D_dash(:,:,2);
            obj.h = rs(2) - rs(1);
            obj.r = (obj.eta + 1)/2 + rs(1)/obj.h;
        end

        function [L0, L1, L2] = stiffnessOp(obj, udof, n)
            % stiffnessOp stiffness operator 
            cn = obj.mat.tensor/obj.mat.tensor(1,2,1,2); % normalized stiffness tensor
            % relevant material matrices: 
            czz = squeeze(cn(1,udof,udof,1));
            crr = squeeze(cn(2,udof,udof,2));
            cpp = squeeze(cn(3,udof,udof,3));
            Crp = squeeze(cn(2,udof,udof,3)) + squeeze(cn(3,udof,udof,2));
            Czp = squeeze(cn(1,udof,udof,3)) + squeeze(cn(3,udof,udof,1));
            Crz = squeeze(cn(2,udof,udof,1)) + squeeze(cn(1,udof,udof,2));
            crp = squeeze(cn(2,udof,udof,3));
            crz = squeeze(cn(2,udof,udof,1));
            % differetiation in curvilinear coordinate system:
            A = [0, 0, 0; 0, 0, -1; 0, 1, 0]; A = squeeze(A(udof, udof));
            B = [0, 0, 0; 0, 1,  0; 0, 0, 1]; B = squeeze(B(udof, udof));
            % differentiation matrices on normalized domain:
            Dr1 = obj.D1; Dr2 = obj.D2; r = obj.r;
            Id = eye(size(Dr1)); % identity matrix for discretization
            % operators:
            L2 = kron(czz, Id); 
            L1 = kron(Czp*A, diag(1./r)) + kron(Crz, Dr1) + 1i*n*kron(Czp, diag(1./r));
            L0 = kron(crr, Dr2) + kron(Crp*A, diag(1./r)*Dr1 - diag(1./r.^2)) ...
                + kron(cpp, diag(1./r)*Dr1) - kron(cpp*B, diag(1./r.^2)) ...
                + 1i*n*(kron(Crp, diag(1./r)*Dr1 - diag(1./r.^2)) + kron(2*cpp*A, diag(1./r.^2))) ...
                + (1i*n).^2*(kron(cpp, diag(1./r.^2)));
        end

        function M = massOp(obj, udof)
            % massOp mass operator 
            rhon = 1*eye(length(udof)); % normalized mass matrix tensor 
            Id = eye(size(obj.D1)); % identity matrix for discretization
            M = kron(rhon, Id); 
        end

        function [B0, B1] = tractionOp(obj, udof, n)
            % tractionOp traction operator 
            cn = obj.mat.tensor/obj.mat.tensor(1,2,1,2); % normalized stiffness tensor
            % relevant material matrices: 
            crr = squeeze(cn(2,udof,udof,2));
            crp = squeeze(cn(2,udof,udof,3));
            crz = squeeze(cn(2,udof,udof,1));
            A = [0, 0, 0; 0, 0, -1; 0, 1, 0]; A = squeeze(A(udof, udof)); % differetiation in curvilinear coordinate system
            % differentiation matrices on normalized domain 
            Dr1 = obj.D1; rInv = diag(1./obj.r); 
            Id = eye(obj.N); % identity matrix for discretization
            % operators:
            B1 = kron(crz, Id([1, obj.N], :)); 
            B0 = kron(crr, Dr1([1, obj.N], :)) + kron(crp*A, rInv([1, obj.N], :)) + 1i*n*kron(crp, rInv([1, obj.N], :)); 
        end

        function [U0] = displacementOp(obj, udof)
            % displacementOp displacement operator (boundaries)
            I = eye(length(udof));
            Id = eye(obj.N);  
            U0 = kron(I, Id([1, obj.N], :));
        end
        
    end
end
