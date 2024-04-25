classdef LayerCylindrical < Layer
% LayerCylindrical - Class to represent one layer of a multi-layered cylinder.
% There is usually no need to use this class explicitly (used internally by 
% Cylinder).
%
% See also Cylinder, Waveguide.
% 
% 2022-2024 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France

    properties 
        PPr     % integral of weighted product matrix of ansatz functions ∫P*P r dr
        PPdr    % integral of weighted product matrix of ansatz functions ∫P*P' r dr 
        PdPdr   % integral of weighted product matrix of ansatz functions ∫P'*P' r dr 
        PPInvr  % integral of weighted product matrix of ansatz functions ∫P*P 1/r dr 
        PPInvrr % integral of weighted product matrix of ansatz functions ∫P*P 1/r^2 dr 
        PPdInvr % integral of weighted product matrix of ansatz functions ∫P*P' 1/r dr 
    end
    properties (Dependent)
        r % alias to the nodal coordinates obj.y
    end

    methods
        function obj = LayerCylindrical(mat, rs, N)
            % LayerCylindrical - Create a LayerCylindrical object.
            [yn, ~] = Layer.nodes(N);  % nodal coordinates and integration weights
            h = rs(end)-rs(1); % layer thickness
            if rs(1)/h < 1e-2  % if inner radius is close to zero
                basis = Layer.basisGaussLegendre(yn); % avoid singularity at r = 0
            else
                basis = Layer.basisGaussLobattoLumped(yn); % faster integration, diagonal mass
            end
            obj = obj@Layer(mat, rs, N, basis);
            % % element matrices specific to cylindrical coordinates:
            basis.r = basis.y(:) + rs(1)/h; % map nodes on [0, 1] -> [r1,r2]/(r2-r1)
            obj.PPr = LayerCylindrical.elemPPr(basis.P, basis.w, basis.r);             % ∫P*P r dr
            obj.PPdr = LayerCylindrical.elemPPdr(basis.P, basis.Pd, basis.w, basis.r); % ∫P*P' r dr 
            obj.PdPdr = LayerCylindrical.elemPdPdr(basis.Pd, basis.w, basis.r);        % ∫P'*P' r dr 
            obj.PPInvr = LayerCylindrical.elemPPInvr(basis.P, basis.w, basis.r);       % ∫P*P 1/r dr 
            obj.PPInvrr = LayerCylindrical.elemPPInvrr(basis.P, basis.w, basis.r);     % ∫P*P 1/r^2 dr 
            obj.PPdInvr = LayerCylindrical.elemPPdInvr(basis.P, basis.Pd, basis.w, basis.r); % ∫P*P' 1/r dr 
        end

        function r = get.r(obj)
            r = obj.y; % just another name
        end

        function [L0, L1, L2] = stiffnessOp(obj, udof, hl, n)
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
            cxpA = cxp*(A + 1i*n*I);
            Acpx = (A + 1i*n*I)*cpx; 
            Acpr = (A + 1i*n*I)*cpr; 
            AcppA = (A + 1i*n*I)*cpp*(A + 1i*n*I); 

            % reduce to desired polarization (udof):
            cxx  = squeeze(cxx(udof,udof));
            cxr  = squeeze(cxr(udof,udof));
            cxpA = squeeze(cxpA(udof,udof));
            Acpx = squeeze(Acpx(udof,udof));
            Acpr = squeeze(Acpr(udof,udof));
            AcppA = squeeze(AcppA(udof,udof));
            
            % element stiffness:
            K2 = kron(cxx, obj.PPr)*hl^2;
            K1 = kron(cxr, obj.PPdr)*hl + kron(cxpA + Acpx, obj.PP)*hl;
            K0 = kron(Acpr, obj.PPd) + kron(AcppA, obj.PPInvr);
            % element flux:
            [G0, G1] = obj.tractionOp(udof, hl, n);  % not yet scaled by hl
            % combine to polynomial of (ik):
            L2 = K2; L1 = K1 + G1; L0 = K0 + G0;
        end

        function M = massOp(obj, udof, hl)
            % massOp - mass operator M
            rhon = eye(length(udof)); % normalized mass matrix (for each dof in u) 
            M = kron(rhon, obj.PPr); % assemble
            M = M*hl^2;
        end

        function [G0, G1] = tractionOp(obj, udof, hl, n)
            % tractionOp - traction operators (flux, used internally)
            cn = obj.mat.c/obj.mat.c(1,2,1,2); % normalized stiffness tensor
            % relevant material matrices: 
            crx = squeeze(cn(2,:,:,1));
            crr = squeeze(cn(2,:,:,2));
            crp = squeeze(cn(2,:,:,3));
            A = [0, 0, 0; 0, 0, -1; 0, 1, 0]; % differetiation in curvilinear coordinate system
            I = eye(size(A));
            % terms due to curvature (to be done before reducing to udof!)
            crpA = crp*(A + 1i*n*I);
            % reduce to desired polarization (udof):
            crx   = squeeze(crx(udof,udof));
            crr   = squeeze(crr(udof,udof));
            crpA  = squeeze(crpA(udof,udof));
            % normalized element flux:
            G1 = kron( crx , -obj.PPdr.' )*hl;
            G0 = kron( crr , -obj.PdPdr.' )  +  kron( crpA , -obj.PPd.' );
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
            cxpA = cxp*(A + 1i*n*I);
            Acpx = (A + 1i*n*I)*cpx; 
            Acpr = (A + 1i*n*I)*cpr; 
            AcppA = (A + 1i*n*I)*cpp*(A + 1i*n*I); 
            crpA = crp*(A + 1i*n*I);

            % define dofs and continuous operator coefficients:
            xy = [1 2];    % Lamb polarization: flexural and longitudinal
            z =  3;        % SH polarization: torsional
            L2 = abs(cxx); % complex values: we only want to know if certain component are zero or not
            L1 = abs(cxr + crx + cxpA + Acpx);
            L0 = abs(AcppA + Acpr + crr + crpA);
            
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
        function ppri = elemPPInvr(P, w, r) 
            % elemPPInvr - integral ∫P*P*1/r dr of basis functions P
            if r(1) <= 10*eps
                error('GEWTOOL:LayerCylindrical',...
                    'Integrating at r=0 is not possible due to the singularity of 1/r. Switch integration scheme.');
            end
            PtimesP = P.*permute(P,[1 3 2]); % size: [integration points, len P, len P]
            PtimesPr = 1./r.*PtimesP;
            ppri = squeeze( sum(w.'.*PtimesPr,1) );
        end
        function ppri2 = elemPPInvrr(P, w, r) 
            % elemPPInvr - integral ∫P*P*1/r^2 dr of basis functions P
            if r(1) <= 10*eps
                error('GEWTOOL:LayerCylindrical',...
                    'Integrating at r=0 is not possible due to the singularity of 1/r. Switch integration scheme.');
            end
            PtimesP = P.*permute(P,[1 3 2]); % size: [integration points, len P, len P]
            PtimesPr2 = 1./(r.^2).*PtimesP;
            ppri2 = squeeze( sum(w.'.*PtimesPr2,1) );
        end
        function ppdri = elemPPdInvr(P, Pd, w, r) 
            % elemPPInvr - integral ∫P*P'*1/r dr of basis functions P
            if r(1) <= 10*eps
                error('GEWTOOL:LayerCylindrical',...
                    'Integrating at r=0 is not possible due to the singularity of 1/r. Switch integration scheme.');
            end
            PtimesPd = P.*permute(Pd,[1 3 2]); % size: [integration points, len P, len P]
            PtimesPdri = 1./r.*PtimesPd;
            ppdri = squeeze( sum(w.'.*PtimesPdri,1) );
        end
        function ppr = elemPPr(P, w, r) 
            % elemPPr - integral ∫P*P*r dr of basis functions P
            PtimesP = P.*permute(P,[1 3 2]); % size: [integration points, len P, len P]
            PtimesPr = r.*PtimesP;
            ppr = squeeze( sum(w.'.*PtimesPr,1) );
        end
        function ppdr = elemPPdr(P, Pd, w, r) 
            % elemPPdr - integral ∫P*P'*r dr of basis functions P
            PtimesPd = P.*permute(Pd,[1 3 2]); % size: [integration points, len P, len P]
            PtimesPdr = r.*PtimesPd;
            ppdr = squeeze( sum(w.'.*PtimesPdr,1) );
        end
        function ppdr = elemPdPdr(Pd, w, r) 
            % elemPdPdr - integral ∫P'*P'*r dr of basis functions P
            PdtimesPd = Pd.*permute(Pd,[1 3 2]); % size: [integration points, len P, len P]
            PdtimesPdr = r.*PdtimesPd;
            ppdr = squeeze( sum(w.'.*PdtimesPdr,1) );
        end
    end
end
