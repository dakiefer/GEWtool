classdef LayerCylindrical < Layer
% LayerCylindrical - Class to represent one layer of a multi-layered cylinder.
% There is usually no need to use this class explicitly (used internally by 
% Cylinder).
%
% See also Cylinder, Waveguide.
% 
% 2022 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France

    properties 
        PPr     % integral of weighted product matrix of ansatz functions ∫P*P r dr
        PPInvr  % integral of weighted product matrix of ansatz functions ∫P*P 1/r dr 
        PPdr    % integral of weighted product matrix of ansatz functions ∫P*P' r dr 
        PdPdr   % integral of weighted product matrix of ansatz functions ∫P'*P' r dr 
    end
    properties (Dependent)
        r % alias to the nodal coordinates obj.y
    end

    methods
        function obj = LayerCylindrical(mat, rs, N)
            % LayerCylindrical - Create a LayerCylindrical object.
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
            % stiffnessOp - stiffness operators L0, L1, L2
            cn = obj.mat.c/obj.mat.c(1,2,1,2); % normalized stiffness tensor
            % relevant material matrices: 
            cxx = squeeze(cn(1,:,:,1));
            cpp = squeeze(cn(3,:,:,3));
            cpr = squeeze(cn(3,:,:,2));
            cxp = squeeze(cn(1,:,:,3));
            cpx = squeeze(cn(3,:,:,1));
            cxr = squeeze(cn(1,:,:,2));
            % differetiation in curvilinear coordinate system:
            A = [0, 0, 0; 0, 0, -1; 0, 1, 0]; % differetiation in curvilinear coordinate system
            I = eye(size(A));

            % include terms due to curvature (to be done before reducing to udof!)
            cxpC = (cxp + cpx)*(1i*n*I + A);
            cprC = cpr*(1i*n*I + A);
            cppC = cpp*(1i*n*I + A)^2;

            % reduce to desired polarization (udof):
            cxx  = squeeze(cxx(udof,udof));
            cpp  = squeeze(cpp(udof,udof));
            cpr  = squeeze(cpr(udof,udof));
            cxp  = squeeze(cxp(udof,udof));
            cpx  = squeeze(cpx(udof,udof));
            cxr  = squeeze(cxr(udof,udof));
            cxpC = squeeze(cxpC(udof,udof));
            cprC = squeeze(cprC(udof,udof));
            cppC = squeeze(cppC(udof,udof));
            
            % element stiffness:
            K2 = kron(cxx, obj.PPr);
            K1 = kron(cxr, obj.PPdr) + kron(cxpC, obj.PP);
            K0 = kron(cpp + cprC, obj.PPd) + kron(cppC - cprC, obj.PPInvr);
            % element flux:
            [G0, G1] = obj.tractionOp(udof, n);
            % combine to polynomial of (ik):
            L2 = K2; L1 = K1 + G1; L0 = K0 + G0;
        end

        function M = massOp(obj, udof)
            % massOp - mass operator M
            rhon = eye(length(udof)); % normalized mass matrix (for each dof in u) 
            M = kron(rhon, obj.PPr); % assemble
        end

        function [G0, G1] = tractionOp(obj, udof, n)
            % tractionOp - traction operators (flux, used internally)
            cn = obj.mat.c/obj.mat.c(1,2,1,2); % normalized stiffness tensor
            % relevant material matrices: 
            crx = squeeze(cn(2,:,:,1));
            crr = squeeze(cn(2,:,:,2));
            crp = squeeze(cn(2,:,:,3));
            A = [0, 0, 0; 0, 0, -1; 0, 1, 0]; % differetiation in curvilinear coordinate system
            I = eye(size(A));
            % terms due to curvature (to be done before reducing to udof!)
            crpC = crp*(1i*n*I + A);
            % reduce to desired polarization (udof):
            crx   = squeeze(crx(udof,udof));
            crr   = squeeze(crr(udof,udof));
            crpC  = squeeze(crpC(udof,udof));
            % normalized element flux:
            G1 = kron( crx , -obj.PPdr.' );
            G0 = kron( crr , -obj.PdPdr.' )  +  kron( crpC , -obj.PPd.' );
        end
        
        function decoupl = decouplesLambvsSH(obj,n)
            % decouplesLambvsSH - Tests whether the xr- and phi-polarized waves decouple.
            if length(obj) > 1
                error('GEWTOOL:decouplesLambvsSH:notimplementedyet', 'Multilayer waveguides do not support this operation at the moment.');
            end

            % stiffnessOp stiffness operator 
            cn = obj.mat.c/obj.mat.c(1,2,1,2); % normalized stiffness tensor

            % relevant material matrices: 
            cxx = squeeze(cn(1,:,:,1));
            crr = squeeze(cn(2,:,:,2));
            cpp = squeeze(cn(3,:,:,3));
            cxr = squeeze(cn(1,:,:,2));
            crx = squeeze(cn(2,:,:,1));
            cxp = squeeze(cn(1,:,:,3));
            cpx = squeeze(cn(3,:,:,1));
            crp = squeeze(cn(2,:,:,3));
            cpr = squeeze(cn(3,:,:,2));

            % differetiation in curvilinear coordinate system:
            A = [0, 0, 0; 0, 0, -1; 0, 1, 0];
            I = eye(size(A));
            cxpC = (cxp + cpx)*(n*I + A); % take n*I real valued for test!
            cprC = cpr*(n*I + A);
            cppC = cpp*(n*I + A)^2;
            crpC = crp*(n*I + A);

            % define dofs and continuous operator coefficients:
            xy = [1 2];   % Lamb polarization: flexural and longitudinal
            z =  3;       % SH polarization: torsional
            L2 = cxx; 
            L1 = cxr + crx + cxpC; 
            L0 = crr + cpp + crpC + cprC + cppC;
            
            % test 
            test = ~(L2(xy,z) | L2(z,xy).' | L1(xy,z) | L1(z,xy).' | L0(xy,z) | L0(z,xy).');
            decoupl = any(test(:));
        end

        %% overload operators: 
        function ret = eq(a, b)
            % eq - Test if the layers a and b are physically the same.
            % True if the layers are of same thickness and material.
            % Usage: 
            % isEq = eq(a, b);
            % isEq = a == b;
            ret = eq@Layer(a, b);
        end
        function ret = ne(a, b)
            ret = ne@Layer(a, b);
        end
        
    end % methods

    methods (Static)
        function ppr = elemPPr(P, w, r) 
            % elemPPr - integral ∫P*P*1/r dr of basis functions P
            PtimesP = P.*permute(P,[1 3 2]); % size: [nodal points, len P, len P]
            PtimesPr = r(:).*PtimesP;
            ppr = squeeze( sum(w.'.*PtimesPr,1) );
        end
        function ppri = elemPPInvr(P, w, r) 
            % elemPPInvr - integral ∫P*P*1/r dr of basis functions P
            PtimesP = P.*permute(P,[1 3 2]); % size: [nodal points, len P, len P]
            PtimesPr = (1./r(:)).*PtimesP;
            ppri = squeeze( sum(w.'.*PtimesPr,1) );
        end
        function ppdr = elemPPdr(P, Pd, w, r) 
            % elemPPdr - integral ∫P*P'*1/r dr of basis functions P
            PtimesPd = P.*permute(Pd,[1 3 2]); % size: [nodal points, len P, len P]
            PtimesPdr = r(:).*PtimesPd;
            ppdr = squeeze( sum(w.'.*PtimesPdr,1) );
        end
        function ppdr = elemPdPdr(Pd, w, r) 
            % elemPdPdr - integral ∫P*P'*1/r dr of basis functions P
            PdtimesPd = Pd.*permute(Pd,[1 3 2]); % size: [nodal points, len P, len P]
            PdtimesPdr = r(:).*PdtimesPd;
            ppdr = squeeze( sum(w.'.*PdtimesPdr,1) );
        end
    end
end
