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
    end
    properties (Dependent)
        r % alias to the nodal coordinates obj.z
    end

    methods
        function obj = LayerCylindrical(mat, rs, N)
            % LayerCylindrical - Create a LayerCylindrical object.
            [rn, ~] = Layer.nodes(N);  % nodal coordinates and integration weights
            h = rs(end)-rs(1); % layer thickness
            if rs(1)/h < 1e-2  % if inner radius is close to zero
                basis = Layer.basisGaussLegendre(rn); % avoid singularity at r = 0
            else
                basis = Layer.basisGaussLobattoLumped(rn); % faster integration, diagonal mass
            end
            obj = obj@Layer(mat, rs, N, basis);
            % % element matrices specific to cylindrical coordinates:
            basis.r = basis.z(:) + rs(1)/h; % map nodes on [0, 1] -> [r1,r2]/(r2-r1)
            obj.PPr = LayerCylindrical.elemPPr(basis.P, basis.w, basis.r);             % ∫P*P r dr
            obj.PPdr = LayerCylindrical.elemPPdr(basis.P, basis.Pd, basis.w, basis.r); % ∫P*P' r dr 
            obj.PdPdr = LayerCylindrical.elemPdPdr(basis.Pd, basis.w, basis.r);        % ∫P'*P' r dr 
            obj.PPInvr = LayerCylindrical.elemPPInvr(basis.P, basis.w, basis.r);       % ∫P*P 1/r dr 
        end

        function r = get.r(obj)
            r = obj.z; % just another name
        end

        function [L0, L1, L2] = stiffnessOp(obj, udof, np, hl, n)
            % stiffnessOp - stiffness operators L0, L1, L2
            
            % relevant material matrices: 
            cn = obj.mat.c/np.c0; % normalized stiffness tensor
            hl = hl/np.h0; % normalize thickness
            cxx = squeeze(cn(1,:,:,1));
            cpp = squeeze(cn(2,:,:,2));
            cpr = squeeze(cn(2,:,:,3));
            cxp = squeeze(cn(1,:,:,2));
            cpx = squeeze(cn(2,:,:,1));
            cxr = squeeze(cn(1,:,:,3));
            crx = squeeze(cn(3,:,:,1)); % boundary flux
            crr = squeeze(cn(3,:,:,3)); % boundary flux
            crp = squeeze(cn(3,:,:,2)); % boundary flux

            % include terms due to curvature (to be done before reducing to udof!)
            A = [0, 0, 0; 0, 0, 1; 0, -1, 0]; % differetiation in curvilinear coordinate system
            I = eye(size(A));
            cxpA = cxp*(A + 1i*n*I);
            Acpx = (A + 1i*n*I)*cpx;
            Acpr = (A + 1i*n*I)*cpr;
            AcppA = (A + 1i*n*I)*cpp*(A + 1i*n*I);
            crpA = crp*(A + 1i*n*I); % flux

            % reduce to desired polarization (udof):
            cxx  = squeeze(cxx(udof,udof));
            cxr  = squeeze(cxr(udof,udof));
            cxpA = squeeze(cxpA(udof,udof));
            Acpx = squeeze(Acpx(udof,udof));
            Acpr = squeeze(Acpr(udof,udof));
            AcppA = squeeze(AcppA(udof,udof));
            crx   = squeeze(crx(udof,udof));  % boundary flux
            crr   = squeeze(crr(udof,udof));  % boundary flux
            crpA  = squeeze(crpA(udof,udof)); % boundary flux
            
            % assemble element stiffness terms:
            K2xx = kron(cxx,   obj.PPr);
            K0pp = kron(AcppA, obj.PPInvr);
            G0rr = kron(crr,  -obj.PdPdr.');
            K1xr = kron(cxr,   obj.PPdr);
            G1xr = kron(crx,  -obj.PPdr.');
            K0pr = kron(Acpr,  obj.PPd);
            G0pr = kron(crpA, -obj.PPd.');
            K1xp = kron(cxpA + Acpx, obj.PP);
            
            % combine to polynomial of (ik):
            L2 = K2xx*hl^2; 
            L1 = (K1xr + G1xr + K1xp)*hl;   % add in this order to preserve hermiticity in numerical presicion!
            L0 = K0pr + G0pr + K0pp + G0rr; % add in this order to preserve hermiticity in numerical presicion!
        end

        function M = massOp(obj, udof, np, hl)
            % massOp - mass operator M
            rhon = obj.mat.rho/np.rho0;
            hl = hl/np.h0; % normalize thickness
            MM = rhon*eye(length(udof)); % normalized mass matrix (for each dof in u) 
            M = kron(MM, obj.PPr)*hl^2; % assemble
        end
        
        function decoupl = decouplesLambvsSH(obj,n)
            % decouplesLambvsSH - Tests whether the x-r- and phi-polarized waves decouple.

            % relevant material matrices: 
            cn = obj.mat.c/obj.mat.c(1,2,1,2); % normalized stiffness tensor
            cxx = squeeze(cn(1,:,:,1));
            crr = squeeze(cn(3,:,:,3));
            cpp = squeeze(cn(2,:,:,2));
            cxr = squeeze(cn(1,:,:,3));
            crx = squeeze(cn(3,:,:,1));
            cxp = squeeze(cn(1,:,:,2));
            cpx = squeeze(cn(2,:,:,1));
            crp = squeeze(cn(3,:,:,2));
            cpr = squeeze(cn(2,:,:,3));

            % include terms due to curvature (to be done before reducing to udof!)
            A = [0, 0, 0; 0, 0, 1; 0, -1, 0]; % differetiation in curvilinear coordinate system
            I = eye(size(A));
            cxpA = cxp*(A + 1i*n*I);
            Acpx = (A + 1i*n*I)*cpx; 
            Acpr = (A + 1i*n*I)*cpr; 
            AcppA = (A + 1i*n*I)*cpp*(A + 1i*n*I); 
            crpA = crp*(A + 1i*n*I);

            % define dofs and continuous operator coefficients:
            lb = Waveguide.udofLamb; % Lamb polarization (axial): flexural and longitudinal
            sh =  Waveguide.udofSH;   % SH polarization (axial): torsional
            iszero = @(cc) all( cc(lb,sh) == 0 ) && all( cc(sh,lb) == 0);

            decoupl = iszero(cxx) && iszero(cxr) && iszero(cxpA + Acpx) && iszero(Acpr) ...
                && iszero(AcppA) && iszero(crx) && iszero(crr) && iszero(crpA);
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
    end
end
