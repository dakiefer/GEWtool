classdef Cylinder < Waveguide
% Cylinder - Represents axially guided waves in cylinders.
% Displacement ansatz: u(x,r,phi,t) = u(r)*exp(i k x + i n phi - i w t)
% 
% Example:
% mat = Material('steel'); % load material data
% r = [5e-3, 6e-3]; % radial coordiantes in m
% N = 20; % discretization (number of nodal points)
% cyl = Cylinder(mat, r, N);
%
% See also Cylinder.Cylinder, Plate, Waveguide.
% 
% 2022-2023 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France

methods 
	function obj = Cylinder(mats, rs, Ns)
        % Cylinder - Create a cylindrical waveguide (axial propagation).
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
        % cyl = Cylinder(mat, r, N);
        %
        % See also: Cylinder, Plate.
		obj = obj@Waveguide(mats, rs, Ns);% converts mats 
        obj.lay = LayerCylindrical.empty; % initialize with correct class
		for ii = 1:length(obj.mat)
			obj.lay(ii) = LayerCylindrical(obj.mat(ii), rs(ii:ii+1), Ns(ii));
		end
    end

    function gew = longitudinal(obj)
        % LONGITUDINAL - Assemble operators for the longitudinal/compressional modes.
        % The longitudional modes (also called compressional) are ux-ur-polarized 
        % and are akin to Lamb waves in a plate. The modes might exist only for
        % circumferential order zero (n=0). Otherwise they are coupled with the
        % torsional waves. They are usually labeled L(0,m).
        % 
        % See also: torsional, flexural, decouplesLambvsSH.
        n = 0;
        gew = obj.Lamb(n);
    end

    function gew = torsional(obj)
        % TORSIONAL - Assemble operators for the torsional modes.
        % The torsional modes are polarized in uPhi and are akin to SH waves in
        % a plate. The modes might exist only for circumferential order zero
        % (n=0). Otherwise they are coupled with the longitudinal waves. They
        % are usually labeled as T(0,m).
        % 
        % See also: longitudinal, flexural, decouplesLambvsSH.
        n = 0;
        gew = obj.sh(n);
    end

    function gew = flexural(obj,n)
        % FLEXURAL - Assemble operators for the flexural modes of circumf. order n.
        % The flexural modes are fully polarized in ux-ur-uPhi and are akin to
        % fully coupled waves in a plate. They are usually labeled F(n,m).
        % 
        % See also: torsional, longitudinal, decouplesLambvsSH.
        gew = obj.fullyCoupled(n);
    end

end % methods

end % class
