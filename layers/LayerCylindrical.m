classdef LayerCylindrical < Layer
    properties (Dependent)
        r
    end

    methods
        function obj = LayerCylindrical(mat, rs, N)
            % LayerCylindrical: constructor
            obj = obj@Layer(mat, rs, N);
        end

        function r = get.r(obj)
            r = obj.y; % just another name
        end

        function [L0, L1, L2] = stiffnessOp(obj, udof, n)
            % stiffnessOp stiffness operator 
            cn = obj.mat.c/obj.mat.c(1,2,1,2); % normalized stiffness tensor
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

        function [B0, B1] = tractionOp(obj, udof, n)
            % tractionOp traction operator 
            cn = obj.mat.c/obj.mat.c(1,2,1,2); % normalized stiffness tensor
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
        
    end
end
