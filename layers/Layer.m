classdef Layer 

    properties (Access = public)
        y   % nodal points in physical domain
        eta % normalized nodal points
        h   % thickness in m
        N   % number of collocation points
        mat % material
        PP  % integral of product matrix of ansatz functions ∫P*Pdy
        PPd % integral of product matrix of ∫P*P'dy
        PdPd % integral of product matrix of ∫P'*P'dy
        D1  % diff matrix on unit domain 
    end

    methods
        function obj = Layer(mat, ylim, N)
            % LayerCylindrical: constructor

            % % element matrices:
            [yi, wi] = Layer.nodes(N);  % nodal coordinates and integration weights
            [P, Pd] = Layer.basis(yi, N);  % polynomial basis
            obj.D1 = collocD(yi);       % differentiation matrix for given basis and nodes
            obj.PP = Layer.elemPP(P, wi);
            obj.PPd = Layer.elemPPd(P, Pd, wi);
            obj.PdPd = Layer.elemPdPd(Pd, wi);

            % % save coordinates and other properties:
            obj.h = ylim(end) - ylim(1); % physical thickness
            obj.y = obj.h*yi + ylim(1);
            obj.eta = yi;
            obj.mat = mat; % material
            obj.N = N; % polynomial order
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
%             [yi, wi] = chebpts(2*N, [0 1]); % integration weights: integrate exactly P*P
            % % For Lagrange polynomials with GLL points:
            [yi, wi] = lglnodes(N-1); yi = flip(yi); wi = wi.';
            wi = wi/2; yi = yi/2 + 1/2; % scale to [0, 1]
        end

        function [P, Pd] = basis(yi, N)
            if nargin < 2
                N = length(yi); % polynomial order is equal to number of quadrature points
            end
            % % For Chebyshev polynomials: 
%             Dy = diffmat(2*N, [0 1]); % differentiation matrix on integration grid yi
%             Psi = chebpoly(0:N-1, [0 1]); % Chebyshev polynomials
%             P = Psi(yi,:); % along 1st dim: samples at yi, along 2nd dim: polynomial order
%             Pd = squeeze(sum(Dy.*shiftdim(P, -1), 2)); % differentiated polynomials
            % % Lagrange polynomials:
            P = eye(length(yi)); % Psi(yi,:);
            Dy = collocD(yi);
            Pd = squeeze(sum(Dy.*shiftdim(P, -1), 2)); % differentiated polynomials;
        end
    end

    methods (Abstract)
        [L0, L1, L2] = stiffnessOp(obj, udof, n)
        [B0, B1] = tractionOp(obj, udof, n)
    end

end % class
