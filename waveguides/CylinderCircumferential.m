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
    end

    function decoupl = decouplesLambvsSH(obj, ~)
        n = 1; % value is irrelevant for circumferential waves
        decoupl = decouplesLambvsSH@Waveguide(obj, n);
    end

end % methods

methods (Static)
    function udof = udofLamb()
        % udofLamb - return the displacement degrees of freedom of Lamb polarization 
        udof = [2 3]; 
    end
    function udof = udofSH()
        % udofSH - return the displacement degrees of freedom of SH polarization 
        udof = 1; 
    end
end

end % class
