classdef LayerCylindrical < Layer
    properties 
        PPr     % integral of weighted product matrix of ansatz functions ∫1/r*P*P dr
        PPr2    % integral of weighted product matrix of ansatz functions ∫1/r^2*P*P dr 
        PPdr    % integral of weighted product matrix of ansatz functions ∫1/r*P*P' dr 
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
            obj.PPr = LayerCylindrical.elemPPr(P, wi, rn);          % ∫1/r*P*P dr
            obj.PPr2 = LayerCylindrical.elemPPr2(P, wi, rn);        % ∫1/r^2*P*P dr 
            obj.PPdr = LayerCylindrical.elemPPdr(P, Pd, wi, rn);    % ∫1/r*P*P' dr 
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
            crp = squeeze(cn(2,udof,udof,3));
            cpr = squeeze(cn(3,udof,udof,2));
            cxp = squeeze(cn(1,udof,udof,3));
            cpx = squeeze(cn(3,udof,udof,1));
            cxr = squeeze(cn(1,udof,udof,2));
            % differetiation in curvilinear coordinate system:
            A = [0, 0, 0; 0, 0, -1; 0, 1, 0]; A = squeeze(A(udof, udof));
            B = [0, 0, 0; 0, 1,  0; 0, 0, 1]; B = squeeze(B(udof, udof));
            I = eye(size(A));
            
            % element stiffness:
            K2 = kron(cxx, obj.PP);
            K1 = kron(cxr, obj.PPd) + kron( (cxp + cpx)*(1i*n*I + A) , obj.PPr);
            k0PPdr = cpp + cpr*(1i*n*I + A); % first term in K0
            k0PPr2 = -cpp*B - (crp + cpr)*A + (1i*n)*(cpp*2*A - (crp + cpr)*I) + (1i*n)^2*cpp; % second term in K0
            K0 = kron(k0PPdr, obj.PPdr) + kron(k0PPr2, obj.PPr2);
            % element flux:
            [G0, G1] = obj.tractionOp(udof, n);
            % combine to polynomial of (ik):
            L2 = K2; L1 = K1 + G1; L0 = K0 + G0;
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
            G1 = kron( crx , -obj.PPd.' );
            G0 = kron( crr , -obj.PdPd.' )  +  kron( crp*(1i*n*I + A) , -obj.PPdr.' );
        end
        
    end % methods

    methods (Static)
        function ppr = elemPPr(P, w, r) 
            % elemPPr: integral ∫P*P*1/r dr of basis functions P
            PtimesP = P.*permute(P,[1 3 2]); % size: [nodal points, len P, len P]
            PtimesPr = (1./r(:)).*PtimesP;
            ppr = squeeze( sum(w.'.*PtimesPr,1) );
        end
        function ppr2 = elemPPr2(P, w, r) 
            % elemPPr2: integral ∫P*P*1/r^2 dr of basis functions P
            PtimesP = P.*permute(P,[1 3 2]); % size: [nodal points, len P, len P]
            PtimesPr = (1./r(:).^2).*PtimesP;
            ppr2 = squeeze( sum(w.'.*PtimesPr,1) );
        end
        function ppdr = elemPPdr(P, Pd, w, r) 
            % elemPPdr: integral ∫P*P'*1/r dr of basis functions P
            PtimesPd = P.*permute(Pd,[1 3 2]); % size: [nodal points, len P, len P]
            PtimesPdr = (1./r(:)).*PtimesPd;
            ppdr = squeeze( sum(w.'.*PtimesPdr,1) );
        end
    end
end
