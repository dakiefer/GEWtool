classdef LayerCylCircumferential < LayerCylindrical
% LayerCylCircumferential - Circumferential waves in a layer of a tube/rod.
% There is usually no need to use this class explicitly (used internally by 
% CylinderCircumferential).
%
% See also CylinderCircumferential, Waveguide.
% 
% 2024 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France

methods
    function obj = LayerCylCircumferential(mat, rs, N)
        % LayerCylCircumferential - Create a LayerCylCircumferential object.
        obj = obj@LayerCylindrical(mat,rs,N);
    end

    function [L0, L1, L2] = stiffnessOp(obj, udof, ~, ~)
        % stiffnessOp - stiffness operators L0, L1, L2
        cn = obj.mat.c/obj.mat.c(1,2,1,2); % normalized stiffness tensor
        % relevant material matrices: 
        cpp = squeeze(cn(3,:,:,3));
        cpr = squeeze(cn(3,:,:,2));

        % include terms due to curvature (to be done before reducing to udof!)
        A = [0, 0, 0; 0, 0, -1; 0, 1, 0]; % differetiation in curvilinear coordinate system
        Acpr = A*cpr; 
        AcppA = A*cpp*A; 
        cppA = cpp*A; 
        Acpp = A*cpp;

        % reduce to desired polarization (udof):
        cpr = squeeze(cpr(udof,udof));
        Acpr = squeeze(Acpr(udof,udof));
        cpp = squeeze(cpp(udof,udof));
        AcppA = squeeze(AcppA(udof,udof));
        cppA = squeeze(cppA(udof,udof));
        Acpp = squeeze(Acpp(udof,udof));
        
        % element stiffness:
        K2 = kron(cpp, obj.PPInvr);
        K1 = kron(cpr, obj.PPd) + kron(cppA + Acpp, obj.PPInvr);
        K0 = kron(Acpr, obj.PPd) + kron(AcppA, obj.PPInvr);
        % element flux:
        [G0, G1] = obj.tractionOp(udof); % not yet scaled by hl
        % combine to polynomial of (in):
        L2 = K2; L1 = K1 + G1; L0 = K0 + G0;
    end

    function M = massOp(obj, udof, hl)
        % massOp - mass operator M
        rhon = eye(length(udof)); % normalized mass matrix (for each dof in u) 
        M = kron(rhon, obj.PPr); % assemble
        M = M*hl^2;
    end

    function [G0, G1] = tractionOp(obj, udof, ~, ~)
        % tractionOp - traction operators (flux, used internally)
        cn = obj.mat.c/obj.mat.c(1,2,1,2); % normalized stiffness tensor
        % relevant material matrices: 
        crr = squeeze(cn(2,:,:,2));
        crp = squeeze(cn(2,:,:,3));
        A = [0, 0, 0; 0, 0, -1; 0, 1, 0]; % differetiation in curvilinear coordinate system
        % terms due to curvature (to be done before reducing to udof!)
        crpA = crp*A;
        % reduce to desired polarization (udof):
        crr   = squeeze(crr(udof,udof));
        crp  = squeeze(crp(udof,udof));
        crpA  = squeeze(crpA(udof,udof));
        % normalized element flux:
        G1 = kron( crp , -obj.PPd.' );
        G0 = kron( crr , -obj.PdPdr.' )  +  kron( crpA , -obj.PPd.' );
    end
    
    function decoupl = decouplesLambvsSH(obj, ~)
        % decouplesLambvsSH - Tests whether the r-phi- and x-polarized waves decouple.

        % relevant material matrices: 
        cn = obj.mat.c/obj.mat.c(1,2,1,2); % normalized stiffness tensor
        crr = squeeze(cn(2,:,:,2));
        cpp = squeeze(cn(3,:,:,3));
        crp = squeeze(cn(2,:,:,3));
        cpr = squeeze(cn(3,:,:,2));

        % include terms due to curvature
        A = [0, 0, 0; 0, 0, -1; 0, 1, 0]; % differetiation in curvilinear coordinate system
        Acpr = A*cpr; 
        AcppA = A*cpp*A; 
        cppA = cpp*A; 
        Acpp = A*cpp;
        crpA = crp*A;

        % define dofs and continuous operator coefficients:
        rp = [2 3];   % circumferential Lamb polarization
        x =  1;       % circumferential SH polarization
        iszero = @(cc) all( cc(rp,x) == 0 ) && all( cc(x,rp) == 0);

        decoupl = iszero(cpp) && iszero(cpr) && iszero(crp) && iszero(crr)...
            && iszero(cppA + Acpp) && iszero(Acpr) && iszero(AcppA) && iszero(crpA);
    end

end % methods

end
