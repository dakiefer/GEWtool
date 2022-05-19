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
            cxx = squeeze(cn(1,udof,udof,1));
            crr = squeeze(cn(2,udof,udof,2));
            cpp = squeeze(cn(3,udof,udof,3));
            Crp = squeeze(cn(2,udof,udof,3)) + squeeze(cn(3,udof,udof,2));
            Cxp = squeeze(cn(1,udof,udof,3)) + squeeze(cn(3,udof,udof,1));
            Crx = squeeze(cn(2,udof,udof,1)) + squeeze(cn(1,udof,udof,2));
            % differetiation in curvilinear coordinate system:
            A = [0, 0, 0; 0, 0, -1; 0, 1, 0]; A = squeeze(A(udof, udof));
            B = [0, 0, 0; 0, 1,  0; 0, 0, 1]; B = squeeze(B(udof, udof));
            % differentiation matrices on normalized domain:
            Dr1 = obj.D1; Dr2 = obj.D2; r = obj.r;
            Id = eye(size(Dr1)); % identity matrix for discretization
            % operators:
            L2 = kron(cxx, Id); 
            L1 = kron(Cxp*A, diag(1./r)) + kron(Crx, Dr1) + 1i*n*kron(Cxp, diag(1./r));
            L0 = kron(crr, Dr2) + kron(Crp*A, diag(1./r)*Dr1 - diag(1./r.^2)) ...
                + kron(cpp, diag(1./r)*Dr1) - kron(cpp*B, diag(1./r.^2)) ...
                + 1i*n*(kron(Crp, diag(1./r)*Dr1 - diag(1./r.^2)) + kron(2*cpp*A, diag(1./r.^2))) ...
                + (1i*n).^2*(kron(cpp, diag(1./r.^2)));
        end

        function [B0, B1] = tractionOp(obj, udof, n)
            % tractionOp traction operator 
            cn = obj.mat.c/obj.mat.c(1,2,1,2); % normalized stiffness tensor
            % relevant material matrices: 
            crx = squeeze(cn(2,udof,udof,1));
            crr = squeeze(cn(2,udof,udof,2));
            crp = squeeze(cn(2,udof,udof,3));
            A = [0, 0, 0; 0, 0, -1; 0, 1, 0]; A = squeeze(A(udof, udof)); % differetiation in curvilinear coordinate system
            % differentiation matrices on normalized domain 
            Dr1 = obj.D1; rInv = diag(1./obj.r); 
            Id = eye(obj.N); % identity matrix for discretization
            % operators:
            B1 = kron(crx, Id([1, obj.N], :)); 
            B0 = kron(crr, Dr1([1, obj.N], :)) + kron(crp*A, rInv([1, obj.N], :)) + 1i*n*kron(crp, rInv([1, obj.N], :)); 
        end
        
        function decoupl = decouplesLambvsSH(obj)
            if length(obj) > 1
                error('GUWTool:decouplesLambvsSH() needs a scalar layer object.');
            end
            % stiffnessOp stiffness operator 
            cn = obj.mat.c/obj.mat.c(1,2,1,2); % normalized stiffness tensor
            % relevant material matrices: 
            cxx = squeeze(cn(1,:,:,1));
            crr = squeeze(cn(2,:,:,2));
            cpp = squeeze(cn(3,:,:,3));
            Crp = squeeze(cn(2,:,:,3)) + squeeze(cn(3,:,:,2));
            Cxp = squeeze(cn(1,:,:,3)) + squeeze(cn(3,:,:,1));
            Crx = squeeze(cn(2,:,:,1)) + squeeze(cn(1,:,:,2));
            % differetiation in curvilinear coordinate system:
            A = [0, 0, 0; 0, 0, -1; 0, 1, 0]; A = squeeze(A(:,:));
            B = [0, 0, 0; 0, 1,  0; 0, 0, 1]; B = squeeze(B(:,:));
            
            % define dofs and continuous operator coefficients:
            xy = [1 2]; % Lamb polarization
            z =  3; % SH polarization 
            L2 = cxx;
            L1 = Cxp*A + Crx + Cxp;
            L0 = crr + Crp*A + 2*cpp + cpp*B + Crp + 2*cpp*A;
            
            % test 
            test = ~(L2(xy,z) | L2(z,xy).' | L1(xy,z) | L1(z,xy).' | L0(xy,z) | L0(z,xy).');
            decoupl = any(test(:));
        end
        
    end
end
