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

    function [L0, L1, L2] = stiffnessOp(obj, udof, np, ~, ~)
        % stiffnessOp - stiffness operators L0, L1, L2
        cn = obj.mat.c/np.c0; % normalized stiffness tensor
        % relevant material matrices: 
        cpp = squeeze(cn(3,:,:,3));
        cpr = squeeze(cn(3,:,:,2));
        crr = squeeze(cn(2,:,:,2)); % boundary flux
        crp = squeeze(cn(2,:,:,3)); % boundary flux

        % include terms due to curvature (to be done before reducing to udof!)
        A = [0, 0, 0; 0, 0, -1; 0, 1, 0]; % differetiation in curvilinear coordinate system
        Acpr = A*cpr; 
        AcppA = A*cpp*A; 
        cppA = cpp*A; 
        Acpp = A*cpp;
        crpA = crp*A;

        % reduce to desired polarization (udof):
        cpr = squeeze(cpr(udof,udof));
        Acpr = squeeze(Acpr(udof,udof));
        cpp = squeeze(cpp(udof,udof));
        AcppA = squeeze(AcppA(udof,udof));
        cppA = squeeze(cppA(udof,udof));
        Acpp = squeeze(Acpp(udof,udof));
        crr   = squeeze(crr(udof,udof));  % boundary flux
        crp  = squeeze(crp(udof,udof));   % boundary flux
        crpA  = squeeze(crpA(udof,udof)); % boundary flux
        
        % assemble element stiffness terms:
        K2pp = kron(cpp,   obj.PPInvr);
        K1pp = kron(cppA + Acpp, obj.PPInvr);
        K0pp = kron(AcppA, obj.PPInvr);
        G0rr = kron(crr,  -obj.PdPdr.');
        K1pr = kron(cpr,   obj.PPd);
        G1pr = kron(crp,  -obj.PPd.');
        K0pr = kron(Acpr,  obj.PPd);
        G0pr = kron(crpA, -obj.PPd.');

        % combine to polynomial of (in):
        L2 = K2pp;
        L1 = K1pr + G1pr + K1pp;        % add in this order to preserve hermiticity in numerical precision!
        L0 = K0pr + G0pr + G0rr + K0pp; % add in this order to preserve hermiticity in numerical precision!
    end

    function decoupl = decouplesPolarization(obj, dof, ~)
        % decouplesPolarization - Tests whether 'dof' decouples from the remaining 
        % degrees of freedom. 
        rem = setdiff(1:3, dof); % remaining degrees of freedom 
        if isempty(rem), decoupl = true; return; end % no need to test

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
        iszero = @(c) all( c(dof,rem) == 0, 'all' ) && all( c(rem,dof) == 0, 'all' );
        decoupl = iszero(cpp) && iszero(cpr) && iszero(crp) && iszero(crr)...
            && iszero(cppA + Acpp) && iszero(Acpr) && iszero(AcppA) && iszero(crpA);
    end

end % methods

end
