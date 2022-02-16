classdef Layer 

    properties (Access = public)
        y   % collocation points on domain of unit length
        eta % collocation points on [-1, 1] TODO: is this being used?
        h   % thickness in m
        N   % number of collocation points
        D1  % diff matrix on unit domain 
        D2  % second order diff matrix
        mat % material 
    end

    methods
        function obj = Layer(mat, ys, N)
            % LayerCylindrical: constructor
            obj.mat = mat;
            obj.N = N;
            [obj.eta, D_dash] = chebdif(obj.N, 2);
            obj.D1 = 2*D_dash(:,:,1);
            obj.D2 = 4*D_dash(:,:,2);
            obj.h = ys(2) - ys(1);
            obj.y = (obj.eta + 1)/2 + ys(1)/obj.h;
        end

        function M = massOp(obj, udof)
            % massOp mass operator 
            rhon = eye(length(udof)); % normalized mass matrix tensor 
            Id = eye(size(obj.D1)); % identity matrix for discretization
            M = kron(rhon, Id); 
        end

        function U0 = displacementOp(obj, udof)
            % displacementOp displacement operator (boundaries)
            I = eye(length(udof));
            Id = eye(obj.N);  
            U0 = kron(I, Id([1, obj.N], :));
        end
    end % methods

    methods (Abstract)
        [L0, L1, L2] = stiffnessOp(obj, udof, n)
        [B0, B1] = tractionOp(obj, udof, n)
    end

end % class
