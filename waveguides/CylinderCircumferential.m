classdef CylinderCircumferential < Cylinder
% CylinderCircumferential - Represents circumferential guided waves in cylinders.
% Displacement ansatz: u(x,phi,r,t) = u(r)*exp(i k phi b - i w t), where: 
% - b: outer radius 
% - phi: angular coordinate 
% - k wavenumber at the outer cylinder's surface
% 
% Example:
% mat = Material('steel'); % load material data
% r = [5e-3, 6e-3]; % radial coordiantes in m
% N = 20; % discretization (number of nodal points)
% cyl = CylinderCircumferential(mat, r, N);
%
% See also Cylinder.Cylinder, Plate, Waveguide.
% 
% 2024 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France

methods 
	function obj = CylinderCircumferential(mats, rs, Ns)
        % CylinderCircumferential - Cylindrical waveguide with circumferential wave propagation.
        % Arguments: 
        % - mats:  materials [1 x Nlay], either of class "Material" or a struct
        %          needs to support mats.rho (scalar) and mats.c (3x3x3x3).
        % - rs:    radial coordinates of interfaces in meter [1 x Nlay+1]
        % - Ns:    discretization order for each layer [1 x Nlay]
        % 
        % Example:
        % mat = Material('steel'); % load material data
        % r = [5e-3, 6e-3]; % radial coordiantes in m
        % N = 20; % discretization (number of nodal points)
        % circ = CylinderCircumferential(mat, r, N);
        %
        % See also: Cylinder, Plate.
		obj = obj@Cylinder(mats, rs, Ns); % converts mats 
        for ii = 1:length(obj.mat)
			obj.lay{ii} = LayerCylCircumferential(obj.mat{ii}, rs(ii:ii+1), obj.geom.N(ii));
        end
        obj.n = 0;  % The wavenumber along the axis of the cylinder is always zero.
    end

    function gew = Lamb(obj)
        % LAMB - Assemble operators for uphi-ur-polarized waves.
        % The modes are akin to Lamb waves in a plate. The modes exist only when
        % the material has appropriate symmetries. Otherwise they are coupled
        % with the sh-type waves that are polarized in ux-direction.
        % 
        % See also: sh, fullyCoupled.
        udof = obj.udofLamb();
        gew = obj.polarization(udof, 0); % vanishing axial wavenumber 
        % renormalization: instead of n, use k = n/b at outer radius b:
        b = gew.geom.zItf(end,end)/gew.np.h0; % normalized outer radius 
        gew.op.L2 = gew.op.L2*b^2;
        gew.op.L1 = gew.op.L1*b;
        gew.family = 'Lamb-like';
    end

    function gew = sh(obj)
        % SH - Assemble operators for ux-polarized waves.
        % Purely polarized in axial direction, akin to SH waves in a plate. The
        % modes exist only when the material exhibits appropriate symmetries.
        % Otherwise they are coupled with the Lamb-polarized waves.
        % 
        % See also: fullyCoupled, Lamb.
        udof = obj.udofSH(); % ux component
        gew = obj.polarization(udof, 0); % vanishing axial wavenumber 
        % renormalization: instead of n, use k = n/b at outer radius b:
        b = gew.geom.zItf(end,end)/gew.np.h0; % normalized outer radius 
        gew.op.L2 = gew.op.L2*b^2;
        gew.op.L1 = gew.op.L1*b;
        gew.family = 'SH-like';
    end

    function gew = fullyCoupled(obj)
        % FULLYCOUPLED - Assemble operators of fully-polarized waves.
        % The the waves are polarized in ux-ur-uPhi and are akin to
        % fully coupled waves in a plate. 
        % 
        % See also: Lamb, sh.
        udof = [1 2 3]; % ux, uphi and ur components
        gew = obj.polarization(udof, 0); % vanishing axial wavenumber 
        % renormalization: instead of n, use k = n/b at outer radius b:
        b = gew.geom.zItf(end,end)/gew.np.h0; % normalized outer radius 
        gew.op.L2 = gew.op.L2*b^2;
        gew.op.L1 = gew.op.L1*b;
        gew.family = 'all';
    end

    function F = displGrad(obj, dat)
        % displGrad - Displacement gradient of the provided field.
        udof = obj.udof;
        F = cell(obj.geom.nLay, 1); % allocate for each layer
        A = Cylinder.AphiDerivative; % differetiation in curvilinear coordinate 
        n = dat.k*obj.geom.zItf(end)/obj.np.h0; % k times normalized outer radius ro
        inIpA = 1i*n.*shiftdim(eye(3),-3) + shiftdim(A,-3); % n = dat.k/r (circumferential)
        inIpA = inIpA(:,:,:,:,udof); % reduce according to polarization
        u = displacement(dat); % extract displacements (piezoelec. plate also has potentials)
        for l = 1:obj.geom.nLay
            sizeF = [dat.Nk, dat.Nw, obj.geom.N(l), 3, 3]; % size of F = grad u 
            Fi = zeros(sizeF); % allocate for displacement gradient F = grad u
            lay = obj.lay{l};
            Dr = 1/lay.h*lay.D1*obj.np.h0; % differentiation matrix, normalized
            r = lay.r(:)/obj.np.h0; 
            if r(1) == 0
                r(1) = r(1) + max(r)*(100*eps); % lazyly avoid singularity
            end
            r = shiftdim(r,-2);           % radial coordinates, normalized
            % iku = 1i*dat.n.*u{l}; % ex.F = ik*u (k = dat.n == 0, circumferential waves)
            uu = permute(u{l}, [1, 2, 5, 3, 4]); % [k,w,*,r,u] additional dimension for mult. with diff. mat.
            dru = sum(shiftdim(Dr, -2).*uu, 4); % er.F = ∂u/∂r, dimension 4 is singleton
            uu = permute(u{l}, [1, 2, 3, 5, 4]); % [k,w,r,*,u] additional dimension * for mult. with (in*I + A)
            Au = sum(1./r.*inIpA.*uu,5); % k = n/r
            % Fi(:,:,:,1,udof) = 0;  % ex.F = 0
            Fi(:,:,:,2,:) = Au;      % assign components ephi.F (always has 3 components thanks to inIpA))
            Fi(:,:,:,3,udof) = dru;  % assign components er.F
            F{l} = Fi; 
        end 
    end

end % methods

methods (Static)
    function udof = udofLamb()
        % udofLamb - return the displacement degrees of freedom of Lamb-like 
        % polarization, i.e., [uphi, ur] displacements.
        udof = [2 3]; 
    end
    function udof = udofSH()
        % udofSH - return the displacement degrees of freedom of SH-like 
        % polarization, i.e., ux displacements
        udof = 1; 
    end
end % methods

end % class
