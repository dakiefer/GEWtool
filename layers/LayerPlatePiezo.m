classdef LayerPlatePiezo < LayerPlate
% LayerPlatePiezo - Class to represent one layer of a multi-layered piezoelectric plate.
% There is usually no need to use this class explicitly (used internally by 
% Plate).
%
% See also Plate, Waveguide.
% 
% 2024 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France
    
    methods
        function obj = LayerPlatePiezo(mat, ys, N)
            % LayerPlate - Create a LayerPlate object.
            if ~isa(mat, 'MaterialPiezoelectric')
                error('GEWTOOL:LayerPlatePiezo', 'The material should be of class MaterialPiezoelectric.');
            end
            obj = obj@LayerPlate(mat, ys, N);
        end

        function [M] = massOp(obj, udof, np, hl, ~)
            % assemble mass for all equations: 
            rhon = obj.mat.rho/np.rho0;
            MM  = blkdiag(rhon*eye(length(udof)), 0);  % mechanics only, singular
            M = kron(MM,obj.PP)*hl/np.h0; % assemble
        end
        
        function [L0, L1, L2] = stiffnessOp(obj, udof, np, hl, ~)
            % stiffnessOp - stiffness operators L0, L1, L2
            cn = obj.mat.c/np.c0; % normalized stiffness tensor
            en = obj.mat.e/np.e0;
            En = obj.mat.epsilon/np.eps0; % np.eps0 is average of epsilon tensor components
            hl = hl/np.h0; % normalize
            % relevant material matrices: 
            cxx = squeeze(cn(1,udof,udof,1));
            cxz = squeeze(cn(1,udof,udof,3));
            czx = squeeze(cn(3,udof,udof,1));
            czz = squeeze(cn(3,udof,udof,3));
            exx = squeeze(en(1,1,udof));
            exz = squeeze(en(1,3,udof));
            ezx = squeeze(en(3,1,udof));
            ezz = squeeze(en(3,3,udof));
            Exx = En(1,1);
            Exz = En(1,3);
            Ezx = En(3,1);
            Ezz = En(3,3);
            % assemble equations:
            % equation 1: balance of linear momentum, i.e., div(c:grad(u) + grad(Phi).e) + w^2 rho u = 0
            % equation 2: charge-free material, i.e.,       div(e:grad(u) - eps.grad(Phi)) = 0
            A  = [cxx,    exx;      % ~ (ik)^2
                  exx.'  -Exx];
            B  = [czx,    exz;      % ~ (ik)^1
                  ezx.'  -Ezx];
            Bt = [cxz,    ezx;      % ~ (ik)^1
                  exz.'  -Exz];
            C  = [czz,    ezz;      % ~ (ik)^0
                  ezz.'  -Ezz];
            
            % assemble element stiffness terms:
            K2 = kron(A,  obj.PP); 
            K1 = kron(Bt,  obj.PPd);
            G1 = kron(B, -obj.PPd.');  % boundary flux
            G0 = kron(C, -obj.PdPd.'); % boundary flux
            % combine to polynomial of (ik):
            L0 = G0/hl;  L1 = K1 + G1;  L2 = K2*hl;

            % inhomogeneous Neumann BCs for electrical field (the V-term in the paper):
            if ~isdiag(obj.PP), error('TODO: only implemented for GLL Lagrange elements.'); end
            nTop = size(L1,1); nBot = size(A,1); 
            L1(nTop,nTop) = -1i*obj.mat.eps0/np.eps0;
            L1(nBot,nBot) = +1i*obj.mat.eps0/np.eps0;
        end
        
        function decoupl = decouplesLambvsSH(obj,~)
            % decouplesLambvsSH - Tests whether the Lamb- and SH-polarized waves decouple. Argument n is optional and does nothing.
            warning('TODO this function could be re-used from the LayerPlate class.')
            xy = [1 2]; % Lamb polarization
            z =  3; % SH polarization 
            c1test = obj.mat.c(xy,xy,z,xy); 
            c2test = obj.mat.c(xy,z,xy,xy);
            decoupl = all(c1test(:) == 0) & all(c2test(:) == 0);
        end

        %% overload operators: 
        function ret = eq(a, b)
            ret = eq@Layer(a, b);
        end
        function ret = ne(a, b)
            ret = ne@Layer(a, b);
        end
        
    end % methods
end % class
