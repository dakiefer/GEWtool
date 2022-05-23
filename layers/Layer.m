classdef Layer 

    properties (Access = public)
        y   % collocation points on domain of unit length
        h   % thickness in m
        N   % number of collocation points
        mat % material
        PP  % integral of product matrix of ansatz functions ∫P*Pdy
        PPd % integral of product matrix of ∫P*P'dy
        PdPd % integral of product matrix of ∫P'*P'dy
    end

    methods
        function obj = Layer(mat, ylim, N)
            % LayerCylindrical: constructor

            % properties:
            obj.mat = mat; % material
            obj.N = N; % polynomial order
            obj.y = chebpts(N, ylim); % computing grid: defines polynomial order
            obj.h = ylim(end) - ylim(1); % physical thickness
            
            % % element matrices:
            [yi, wi] = Layer.nodes(N);  % nodal coordinates and integration weights
            [P, Pd] = Layer.basis(yi);  % polynomial basis
            obj.PP = Layer.elemPP(P, wi);
            obj.PPd = Layer.elemPPd(P, Pd, wi);
            obj.PdPd = Layer.elemPdPd(Pd, wi);
        end

        function M = massOp(obj, udof)
            % massOp mass operator 
            rhon = eye(length(udof)); % normalized mass matrix (for each dof in u) 
            me = obj.PP; % element mass
            M = kron(rhon,me); % assemble
        end

%         function U0 = displacementOp(obj, udof)
%             % displacementOp displacement operator (boundaries)
%             I = eye(length(udof));
%             Id = eye(obj.N);  
%             U0 = kron(I, Id([1, obj.N], :));
%         end
    end % methods

    methods (Static)
        function me = elemPP(P, w) 
            % elemPP: integral of product matrix of ansatz functions ∫P*Pdy (element mass)
            PtimesP = P.*permute(P,[1 3 2]);
            me = squeeze( sum(w.'.*PtimesP,1) );
        end
        function le1 = elemPPd(P, Pd, w) 
            % elemPPd: integral of product matrix of ∫P*P'dy (element stiffness and flux)
            PtimesPd = P.*permute(Pd,[1 3 2]);
            le1 = squeeze( sum(w.'.*PtimesPd,1) );
        end
        function g0 = elemPdPd(Pd, w) 
            % elemPdPd: integral of product matrix of ∫P'*P'dy (element flux)
            PdtimesPd = Pd.*permute(Pd,[1 3 2]);
            g0 = squeeze( sum(w.'.*PdtimesPd,1) );
        end

        function [yi, wi] = nodes(N)
            % % For Chebyshev polynomials on Chebyshev points: 
            % [yi, wi] = chebpts(2*N, [0 1]); % integration weights: integrate exactly P*P

            % % For Lagrange polynomials with GLL points:
            [yi, wi] = lobpts(N, [-1, 1]); % does only work on dom = [-1 1]!!!!
            wi = wi/2; yi = yi/2 + 1/2; % scale to [0, 1]
        end

        function [P, Pd] = basis(yi)
            % % For Chebyshev polynomials: 
%             Dy = diffmat(2*N, [0 1]); % differentiation matrix on integration grid yi
%             Psi = chebpoly(0:N-1, [0 1]); % Chebyshev polynomials
%             P = Psi(yi,:); % along 1st dim: samples at yi, along 2nd dim: polynomial order
%             Pd = squeeze(sum(Dy.*shiftdim(P, -1), 2)); % differentiated polynomials

            % % For Lagrange polynomials
            Psi = chebfun.lagrange(yi);
            Psid = diff(Psi);
            P = eye(length(yi)); % Psi(yi,:);
            Pd = Psid(yi,:);
        end

    end

    methods (Abstract)
        [L0, L1, L2] = stiffnessOp(obj, udof, n)
        [B0, B1] = tractionOp(obj, udof, n)
    end

end % class
