classdef Layer 
% Layer - Class to represent one layer of a multi-layered waveguide.
% There is usually no need to use this class explicitly (used internally by 
% Waveguide).
%
% See also Plate, Cylinder, Waveguide.
% 
% 2022-2024 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France

    properties (Access = public)
        z   % nodal points in physical domain
        eta % normalized nodal points
        h   % thickness in m
        N   % number of collocation points
        mat % material
        PP  % integral of product matrix of ansatz functions ∫P*Pdz
        PPd % integral of product matrix of ∫P*P'dz
        PdPd % integral of product matrix of ∫P'*P'dz
        D1  % diff matrix on unit domain 
        w   % integration coeffs on unit domain
    end

    methods
        function obj = Layer(mat, zlim, N, basis)
            % Layer - Create a Layer object.
            % Arguments:
            % - mat:   [1 x 1] material of class "Material" or struct
            %          with mat.rho [1 x 1] and mat.c [3 x 3 x 3 x 3].
            % - zlim:  [1 x 2] coordinates of lower and upper interface.
            % - N:     discretization order [1 x 1] (number of nodal points)
            % - b:     (optional) polynomial basis 
    
            [zn, wn] = Layer.nodes(N);  % Gauss-Lobatto element nodes 
            if nargin < 4 
                basis = Layer.basisGaussLobattoLumped(zn); % default polynomial basis
            end
            obj.D1 = collocD(zn); % save differentiation matrix for later post-processing
            % % element matrices:
            obj.PP = Layer.elemPP(basis.P, basis.w);
            obj.PPd = Layer.elemPPd(basis.P, basis.Pd, basis.w);
            obj.PdPd = Layer.elemPdPd(basis.Pd, basis.w);
            % % save coordinates and other properties:
            obj.h = zlim(end) - zlim(1); % physical thickness
            obj.z = obj.h*zn + zlim(1);
            obj.w = wn; % integration weights on unit domain (post-processing)
            obj.eta = zn;
            obj.mat = mat; % material
            obj.N = N; % number of nodes
        end

        function M = massOp(obj, udof, np, hl)
            % massOp mass operator 
            rhon = obj.mat.rho/np.rho0; 
            MM = rhon*eye(length(udof)); % normalized mass matrix (for each dof in u) 
            M = kron(MM, obj.PP)*hl/np.h0; % assemble
        end
    end % methods

    methods (Static)
        function me = elemPP(P, w) 
            % elemPP - integral ∫P*Pdz of basis functions P (element mass)
            PtimesP = P.*permute(P,[1 3 2]);
            me = squeeze( sum(w.'.*PtimesP,1) );
        end

        function le1 = elemPPd(P, Pd, w) 
            % elemPPd - integral ∫P*P'dz of basis functions P (element stiffness and flux)
            PtimesPd = P.*permute(Pd,[1 3 2]);
            le1 = squeeze( sum(w.'.*PtimesPd,1) );
        end

        function g0 = elemPdPd(Pd, w) 
            % elemPdPd - integral ∫P'*P'dz of basis functions P (element flux)
            PdtimesPd = Pd.*permute(Pd,[1 3 2]);
            g0 = squeeze( sum(w.'.*PdtimesPd,1) );
        end

        function [zn, wn] = nodes(N)
            % nodes - nodal points and integration weights.
            % Gauss-Lobatto points zn on unit domain [0, 1] and the
            % corresponding integration weights.
            [zn, wn] = lglnodes(N-1); zn = flip(zn); wn = wn.';
            wn = wn/2; zn = zn/2 + 1/2; % scale to [0, 1]
        end

        function basis = basisGaussLobattoLumped(zn)
            % basisGaussLobattoLumped - sampled ansatz functions and their derivatives.
            % Returns the N ansatz functions P(z) and their derivatives sampled at
            % the integration nodes zi. Here we use the N element nodes zn € [0,1] as
            % integration nodes zi. This is an inaccurate integration scheme for the
            % polynomials P(z)*P(z) of order 2(N-1) that appear in the mass matrix.
            % Gauss-Lobatto quadrature is accurate up to polynomial order 2*Ni-3,
            % where Ni are the number of integration nodes. The advantage of using
            % Ni = N is that we obtain a diagnoal mass matrix ∫P*Pdz.
            N = length(zn); % number of nodes, polynomial order is N-1
            % % Lagrange polynomials as basis:
            Pn = eye(N); % P(zn,i) ith-Lagrange polynomial Pi(z) sampled at nodes zn
            Dz = collocD(zn); % differentiation matrix
            Pdn = squeeze(sum(Dz.*shiftdim(Pn, -1), 2)); % differentiated Pi sampled at zn
            % % Use the element nodes zn as integration points (zi = zn):
            [zi, wi] = lglnodes(N-1);  % integration nodes and weights
            zi = flip(zi)/2 + 1/2; wi = wi.'/2; % scale to [0, 1]
            basis.P = Pn;   basis.Pd = Pdn;    basis.w = wi;     basis.z = zi;
        end

        function basis = basisGaussLegendre(zn)
            % basisGaussLegendre - sampled ansatz functions and their derivatives. 
            % Returns the N ansatz functions P(z) and their derivatives sampled
            % at the integration nodes zi. Here we use Ni = N Gauss-Legendre
            % integration nodes zi € (0,1), which results in an exact
            % integration scheme. Gauss-Legendre quadrature is exact up to
            % polynomial order 2*Ni - 1. The mass matrix ∫P*Pdz is no longer
            % diagonal. The integration nodes zi exclude the element border,
            % which avoids the singularity at r = 0 in cylindrical coordinates.
            N = length(zn);  % number of nodes, polynomial order is N-1
            Ni = N;  % number of quadrature points
            % % Lagrange polynomials as basis:
            Pn = eye(N); % P(zn,i) ith-Lagrange polynomial Pi(z) sampled at nodes zn
            Dz = collocD(zn); % differentiation matrix
            Pdn = squeeze(sum(Dz.*shiftdim(Pn, -1), 2));  % differentiated Pi sampled at zn
            % % interpolate to integration points:
            [zi, wi] = lgwt(Ni,0,1); % Gauss-Legendre points and weights
            zi = flip(zi); % integration points 
            wi = wi.'; % integration weights 
            Pi = zeros(Ni,N);  % P(zi,i) ith-Lagrange polynomial Pi(z) sampled at nodes zi
            Pdi = zeros(Ni,N); % similar but for Pi'(z) 
            for j = 1:N 
                Pi(:,j)  = barylag([zn,Pn(:,j)],zi);  % interpolate Pj(z) to zi
                Pdi(:,j) = barylag([zn,Pdn(:,j)],zi); % interpolate Pj'(z) to zi
            end
            basis.P = Pi;   basis.Pd = Pdi;    basis.w = wi;    basis.z = zi;
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
        [L0, L1, L2] = stiffnessOp(obj, udof, np, hl, n)
    end

end % class
