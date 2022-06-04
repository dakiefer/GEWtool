classdef LayerCylindrical < Layer
    properties 
        PPr     % integral of weighted product matrix of ansatz functions ∫P*P r dr
        PPInvr  % integral of weighted product matrix of ansatz functions ∫P*P 1/r dr 
        PPdr    % integral of weighted product matrix of ansatz functions ∫P*P' r dr 
        PdPdr   % integral of weighted product matrix of ansatz functions ∫P'*P' r dr 
    end
    properties (Dependent)
        r
    end

    methods
        function obj = LayerCylindrical(mat, rs, N)
            % LayerCylindrical: constructor
            obj = obj@Layer(mat, rs, N);

            % % element matrices specific to cylindrical coordinates:
            [yi, wi] = Layer.nodes(N);  % nodal coordinates and integration weights
            [P, Pd] = Layer.basis(yi);  % polynomial basis
            rn = (obj.r)/obj.h; % radial coordinates normalized to thickness
            obj.PPr = LayerCylindrical.elemPPr(P, wi, rn);          % ∫P*P r dr
            obj.PPInvr = LayerCylindrical.elemPPInvr(P, wi, rn);    % ∫P*P 1/r dr 
            obj.PPdr = LayerCylindrical.elemPPdr(P, Pd, wi, rn);    % ∫P*P' r dr 
            obj.PdPdr = LayerCylindrical.elemPdPdr(Pd, wi, rn);     % ∫P'*P' r dr 
        end

        function r = get.r(obj)
            r = obj.y; % just another name
        end

        function [L0, L1, L2] = stiffnessOp(obj, udof, n)
            % stiffnessOp stiffness operator 
            cn = obj.mat.c/obj.mat.c(1,2,1,2); % normalized stiffness tensor
            % relevant material matrices: 
            cxx = squeeze(cn(1,udof,udof,1));
            cpp = squeeze(cn(3,udof,udof,3));
            cpr = squeeze(cn(3,udof,udof,2));
            cxp = squeeze(cn(1,udof,udof,3));
            cpx = squeeze(cn(3,udof,udof,1));
            cxr = squeeze(cn(1,udof,udof,2));
            % differetiation in curvilinear coordinate system:
            A = [0, 0, 0; 0, 0, -1; 0, 1, 0]; A = squeeze(A(udof, udof));
            I = eye(size(A));
            
            % element stiffness:
            K2 = kron(cxx, obj.PPr);
            K1 = kron(cxr, obj.PPdr) + kron( (cxp + cpx)*(1i*n*I + A) , obj.PP);
            k0PPd = cpp + cpr*(1i*n*I + A); % first term in K0
            k0PPInvr = cpp*(1i*n*I + A)^2 - cpr*(1i*n*I + A); % second term in K0
            K0 = kron(k0PPd, obj.PPd) + kron(k0PPInvr, obj.PPInvr);
            % element flux:
            [G0, G1] = obj.tractionOp(udof, n);
            % combine to polynomial of (ik):
            L2 = K2; L1 = K1 + G1; L0 = K0 + G0;
        end

        function M = massOp(obj, udof)
            % massOp mass operator 
            rhon = eye(length(udof)); % normalized mass matrix (for each dof in u) 
            M = kron(rhon, obj.PPr); % assemble
        end

        function [G0, G1] = tractionOp(obj, udof, n)
            % tractionOp traction operator 
            cn = obj.mat.c/obj.mat.c(1,2,1,2); % normalized stiffness tensor
            % relevant material matrices: 
            crx = squeeze(cn(2,udof,udof,1));
            crr = squeeze(cn(2,udof,udof,2));
            crp = squeeze(cn(2,udof,udof,3));
            A = [0, 0, 0; 0, 0, -1; 0, 1, 0]; A = squeeze(A(udof, udof)); % differetiation in curvilinear coordinate system
            I = eye(size(A));
            % normalized element flux:
            G1 = kron( crx , -obj.PPdr.' );
            G0 = kron( crr , -obj.PdPdr.' )  +  kron( crp*(1i*n*I + A) , -obj.PPd.' );
        end
        
        function decoupl = decouplesLambvsSH(obj)
            if length(obj) > 1
                error('GEWTOOL:decouplesLambvsSH:notimplementedyet', 'Multilayer waveguides do not support this operation at the moment.');
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
            A = [0, 0, 0; 0, 0, -1; 0, 1, 0];
            
            % define dofs and continuous operator coefficients:
            xy = [1 2]; % Lamb polarization
            z =  3; % SH polarization 
            L2 = cxx;
            L1 = Cxp*A + Crx + Cxp;
            L0 = crr + Crp*A + 2*cpp - cpp*A*A + Crp + 2*cpp*A;
            
            % test 
            test = ~(L2(xy,z) | L2(z,xy).' | L1(xy,z) | L1(z,xy).' | L0(xy,z) | L0(z,xy).');
            decoupl = any(test(:));
        end
    end % methods

    methods (Static)
        function ppr = elemPPr(P, w, r) 
            % elemPPr: integral ∫P*P*1/r dr of basis functions P
            PtimesP = P.*permute(P,[1 3 2]); % size: [nodal points, len P, len P]
            PtimesPr = r(:).*PtimesP;
            ppr = squeeze( sum(w.'.*PtimesPr,1) );
        end
        function ppri = elemPPInvr(P, w, r) 
            % elemPPr: integral ∫P*P*1/r dr of basis functions P
            PtimesP = P.*permute(P,[1 3 2]); % size: [nodal points, len P, len P]
            PtimesPr = (1./r(:)).*PtimesP;
            ppri = squeeze( sum(w.'.*PtimesPr,1) );
        end
        function ppdr = elemPPdr(P, Pd, w, r) 
            % elemPPdr: integral ∫P*P'*1/r dr of basis functions P
            PtimesPd = P.*permute(Pd,[1 3 2]); % size: [nodal points, len P, len P]
            PtimesPdr = r(:).*PtimesPd;
            ppdr = squeeze( sum(w.'.*PtimesPdr,1) );
        end
        function ppdr = elemPdPdr(Pd, w, r) 
            % elemPPdr: integral ∫P*P'*1/r dr of basis functions P
            PdtimesPd = Pd.*permute(Pd,[1 3 2]); % size: [nodal points, len P, len P]
            PdtimesPdr = r(:).*PdtimesPd;
            ppdr = squeeze( sum(w.'.*PdtimesPdr,1) );
        end
    end
end
