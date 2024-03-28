classdef Layer 
% Layer - Class to represent one layer of a multi-layered waveguide.
% There is usually no need to use this class explicitly (used internally by 
% Waveguide).
%
% See also Plate, Cylinder, Waveguide.
% 
% 2022-2024 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France

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
        function obj = Layer(mat, ylim, N, basis)
            % Layer - Create a Layer object.
            % Arguments:
            % - mat:   [1 x 1] material of class "Material" or struct
            %          with mat.rho [1 x 1] and mat.c [3 x 3 x 3 x 3].
            % - ylim:  [1 x 2] coordinates of lower and upper interface.
            % - N:     discretization order [1 x 1] (number of nodal points)
            % - b:     (optional) polynomial basis 
    
            [yn, wn] = Layer.nodes(N);  % Gauss-Lobatto element nodes 
            if nargin < 4 
                basis = Layer.basisGaussLobattoLumped(yn); % default polynomial basis
            end
            obj.D1 = collocD(yn); % save differentiation matrix for later post-processing
            % % element matrices:
            obj.PP = Layer.elemPP(basis.P, basis.w);
            obj.PPd = Layer.elemPPd(basis.P, basis.Pd, basis.w);
            obj.PdPd = Layer.elemPdPd(basis.Pd, basis.w);
            % % save coordinates and other properties:
            obj.h = ylim(end) - ylim(1); % physical thickness
            obj.y = obj.h*yn + ylim(1);
            obj.w = wn; % integration weights on unit domain (post-processing)
            obj.eta = yn;
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

        function [yn, wn] = nodes(N)
            % nodes - nodal points and integration weights.
            % Gauss-Lobatto points yn on unit domain [0, 1] and the
            % corresponding integration weights.
            [yn, wn] = lglnodes(N-1); yn = flip(yn); wn = wn.';
            wn = wn/2; yn = yn/2 + 1/2; % scale to [0, 1]
        end

        function basis = basisGaussLobattoLumped(yn)
            % basisGaussLobattoLumped - sampled ansatz functions and their derivatives.
            % Returns the N ansatz functions P(y) and their derivatives sampled at
            % the integration nodes yi. Here we use the N element nodes yn € [0,1] as
            % integration nodes yi. This is an inaccurate integration scheme for the
            % polynomials P(y)*P(y) of order (N-1)^2 that appear in the mass matrix.
            % Gauss-Lobatto quadrature is accurate up to polynomial order 2*Ni-3,
            % where Ni are the number of integration nodes. The advantage of using
            % Ni = N is that we obtain a diagnoal mass matrix ∫P*Pdy.
            N = length(yn); % number of nodes, polynomial order is N-1
            % % Lagrange polynomials as basis:
            Pn = eye(N); % P(yn,i) ith-Lagrange polynomial Pi(y) sampled at nodes yn
            Dy = collocD(yn); % differentiation matrix
            Pdn = squeeze(sum(Dy.*shiftdim(Pn, -1), 2)); % differentiated Pi sampled at yn
            % % Use the element nodes yn as integration points (yi = yn):
            [yi, wi] = lglnodes(N-1);  % integration nodes and weights
            yi = flip(yi)/2 + 1/2; wi = wi.'/2; % scale to [0, 1]
            basis.P = Pn;   basis.Pd = Pdn;    basis.w = wi;     basis.y = yi;
        end

        function basis = basisGaussLegendre(yn)
            % basisGaussLegendre - sampled ansatz functions and their derivatives.
            % Returns the N ansatz functions P(y) and their derivatives sampled at
            % the integration nodes yi. Here we use Ni >= N^2/2-N+1 Gauss-Legendre
            % integration nodes yi € (0,1). Ni is chosen such as to obtain an exact
            % integration scheme. Gauss-Legendre quadrature is exact up to
            % polynomial order 2*Ni - 3. The mass matrix ∫P*Pdy is no longer
            % diagonal. The integration nodes yi exclude the element border, which
            % avoids the singularity at r = 0 in cylindrical coordinates. 
            N = length(yn);  % number of nodes, polynomial order is N-1
            Ni = ceil(N^2/2)-N+1;  % number of quadrature points
            % % Lagrange polynomials as basis:
            Pn = eye(N); % P(yn,i) ith-Lagrange polynomial Pi(y) sampled at nodes yn
            Dy = collocD(yn); % differentiation matrix
            Pdn = squeeze(sum(Dy.*shiftdim(Pn, -1), 2));  % differentiated Pi sampled at yn
            % % interpolate to integration points:
            [yi, wi] = lgwt(Ni,0,1); % Gauss-Legendre points and weights
            yi = flip(yi); % integration points 
            wi = wi.'; % integration weights 
            Pi = zeros(Ni,N);  % P(yi,i) ith-Lagrange polynomial Pi(y) sampled at nodes yi
            Pdi = zeros(Ni,N); % similar but for Pi'(y) 
            for j = 1:N 
                Pi(:,j)  = barylag([yn,Pn(:,j)],yi);  % interpolate Pj(y) to yi
                Pdi(:,j) = barylag([yn,Pdn(:,j)],yi); % interpolate Pj'(y) to yi
            end
            basis.P = Pi;   basis.Pd = Pdi;    basis.w = wi;    basis.y = yi;
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
