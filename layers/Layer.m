classdef Layer 
% Layer - Class to represent one layer of a multi-layered waveguide.
% There is usually no need to use this class explicitly (used internally by 
% Waveguide).
%
% See also Plate, Cylinder, Waveguide.
% 
% 2022 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France

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
        w   % integration coeffs on unit domain
    end

    methods
        function obj = Layer(mat, ylim, N)
            % Layer - Create a Layer object.
            % Arguments:
            % - mat:   [1 x 1] material of class "Material" or struct
            %          with mat.rho [1 x 1] and mat.c [3 x 3 x 3 x 3].
            % - ylim:  [1 x 2] coordinates of lower and upper interface.
            % - N:     discretization order [1 x 1] (number of nodal points)

            % % element matrices:
            [yi, wi] = Layer.nodes(N);  % nodal coordinates and integration weights
            [P, Pd] = Layer.basis(yi, N);  % polynomial basis
            obj.D1 = collocD(yi);       % differentiation matrix for given basis and nodes (post-processing)
            obj.PP = Layer.elemPP(P, wi);
            obj.PPd = Layer.elemPPd(P, Pd, wi);
            obj.PdPd = Layer.elemPdPd(Pd, wi);
            % % save coordinates and other properties:
            obj.h = ylim(end) - ylim(1); % physical thickness
            obj.y = obj.h*yi + ylim(1);
            obj.w = wi; % integration weights on unit domain (post-processing)
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
    end % methods

    methods (Static)
        function me = elemPP(P, w) 
            % elemPP - integral ∫P*Pdy of basis functions P (element mass)
            PtimesP = P.*permute(P,[1 3 2]);
            me = squeeze( sum(w.'.*PtimesP,1) );
        end

        function le1 = elemPPd(P, Pd, w) 
            % elemPPd - integral ∫P*P'dy of basis functions P (element stiffness and flux)
            PtimesPd = P.*permute(Pd,[1 3 2]);
            le1 = squeeze( sum(w.'.*PtimesPd,1) );
        end

        function g0 = elemPdPd(Pd, w) 
            % elemPdPd - integral ∫P'*P'dy of basis functions P (element flux)
            PdtimesPd = Pd.*permute(Pd,[1 3 2]);
            g0 = squeeze( sum(w.'.*PdtimesPd,1) );
        end

        function [yi, wi] = nodes(N)
            % nodes - nodal coordinates and integration weights.

            % % For Chebyshev polynomials on Chebyshev points: 
%             [yi, wi] = chebpts(2*N, [0 1]); % integration weights: integrate exactly P*P
            % % For Lagrange polynomials with GLL points:
            [yi, wi] = lglnodes(N-1); yi = flip(yi); wi = wi.';
            wi = wi/2; yi = yi/2 + 1/2; % scale to [0, 1]
        end

        function [P, Pd] = basis(yi, N)
            % basis - ansatz functions and their derivatives sampled at nodes.
            if nargin < 2
                N = length(yi); % polynomial order is equal to number of quadrature points
            end
            % % For Chebyshev polynomials: 
%             Dy = diffmat(N, [0 1]); % differentiation matrix on integration grid yi
%             Psi = chebpoly(0:N-1, [0 1]); % Chebyshev polynomials
%             P = Psi(yi,:); % along 1st dim: samples at yi, along 2nd dim: polynomial order
%             Pd = squeeze(sum(Dy.*shiftdim(P, -1), 2)); % differentiated polynomials
            % % Lagrange polynomials:
            P = eye(length(yi)); % Psi(yi,:);
            Dy = collocD(yi);
            Pd = squeeze(sum(Dy.*shiftdim(P, -1), 2)); % differentiated polynomials;
        end


        %% overload operators: 
        function ret = eq(a, b)
            % eq - Test if the layers a and b are physically the same.
            % True if the layers are of same thickness and material.
            % Usage: 
            % eq(a, b);
            % a == b;
            ret = a.mat == b.mat && (a.h-b.h)/(a.h+b.h)*2 < 1e4*eps;
        end
        function ret = ne(a, b)
            % ne - Test if the layers a and b are physically different.
            % True if the layers are of different thickness or material.
            % Usage: 
            % ne(a, b);
            % a ~= b;
            ret = ~eq(a, b);
        end

    end

    methods (Abstract)
        [L0, L1, L2] = stiffnessOp(obj, udof, n)
        [B0, B1] = tractionOp(obj, udof, n)
    end

end % class
