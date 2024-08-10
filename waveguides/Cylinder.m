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
% 2022-2024 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France

properties (Access = public)
    n = NaN      % Circumferential wavenumber (entire number)
end

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
        if isscalar(rs) || ( isvector(rs) && length(mats) ~= length(rs)-1 )
            error('GEWTOOL:Cylinder', 'Second argument should be a 2-vector containing the inner radius and outer radius as well as interfaces between the layers, all in ascending order.');
        end
		obj = obj@Waveguide(mats, rs, Ns);% converts mats 
        obj.lay = LayerCylindrical.empty; % initialize with correct class
		for ii = 1:length(obj.mat)
			obj.lay(ii) = LayerCylindrical(obj.mat{ii}, rs(ii:ii+1), obj.geom.N(ii));
		end
    end

    function gew = longitudinal(obj)
        % LONGITUDINAL - Assemble operators for the longitudinal/compressional
        % modes. The longitudional modes (also called compressional) are
        % ux-ur-polarized and are akin to Lamb waves in a plate. The modes exist
        % only for circumferential order zero (n=0) when the material has
        % appropriate symmetries. Otherwise they are coupled with the torsional
        % waves. They are usually labeled L(0,m).
        % 
        % See also: torsional, flexural, decouplesLambvsSH.
        n = 0;
        gew = obj.Lamb(n);
    end

    function gew = torsional(obj)
        % TORSIONAL - Assemble operators for the torsional modes. The torsional
        % modes are polarized in uPhi and are akin to SH waves in a plate. The
        % modes exist only for circumferential order zero (n=0) when the
        % material exhibits appropriate symmetries. Otherwise they are coupled
        % with the longitudinal waves. They are usually labeled as T(0,m).
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

    function decoupl = decouplesLambvsSH(obj,n)
        % decouplesLambvsSH - Test decoupling of longitudinal and torsional waves.
        % 
        % Tests whether the Lamb-like (ux-ur displacement, "longitudinal") and
        % SH-like (uphi, "torsional") polarized waves decouple.
        %
        % Arguments:
        % - obj:  (Cylinder) guided wave description
        % - n:    (integer) circumferential wavenumber. If the wave family was
        %         already chosen, this argument is not nessary and the value saved in
        %         obj.n will be used. The explicitly provided argument takes
        %         presedence.
        % 
        % See also: longitudinal, torsional, flexural.
        if nargin < 2 && ~isnan(obj.n)
            n = obj.n; 
        elseif nargin < 2 && isnan(obj.n)
            error('GEWTOOL:Cylinder:decouplesLambvsSH','You must specify the circumferential order n.');
        end
        decoupl = decouplesLambvsSH@Waveguide(obj,n);
    end

    function F = displGrad(obj, dat)
        % displGrad - Displacement gradient of the provided field "dat.u".
        udof = obj.udof;
        F = cell(obj.geom.nLay, 1); % allocate for each layer
        inIpA = 1i*obj.n*eye(3) + [0, 0, 0; 0, 0, -1; 0, 1, 0]; % differetiation in curvilinear coordinate system
        inIpA = inIpA(:,udof); % reduce according to polarization
        for l = 1:obj.geom.nLay
            sizeF = [size(dat.w), obj.geom.N(l), 3, 3]; % size of F = grad u 
            Fi = zeros(sizeF); % allocate for displacement gradient F = grad u
            lay = obj.lay(l);
            Dr = 1/lay.h*lay.D1; % differentiation matrix
            r = lay.r(:); 
            if r(1) == 0
                r(1) = r(1) + max(r)*(100*eps); % lazyly avoid singularity
            end
            r = shiftdim(r,-2);           % radial coordinates in SI units
            iku = 1i*dat.k.*dat.u{l}; % ex.F = ik*u
            uu = permute(dat.u{l}, [1, 2, 5, 3, 4]); % additional dimension for mult. with diff. mat.
            dru = sum(shiftdim(Dr, -2).*uu, 4); % ey.F = ∂u/∂y, dimension 4 is singleton
            uu = permute(dat.u{l}, [1, 2, 3, 5, 4]); % additional dimension for mult. with (in*I + A)
            Ad = shiftdim(inIpA, -3);
            Au = sum(1./r.*Ad.*uu,5);
            Fi(:,:,:,1,udof) = iku;  % assign components ex.F
            Fi(:,:,:,2,udof) = dru;  % assign components er.F
            Fi(:,:,:,3,:) = Au;   % assign components ephi.F
            F{l} = Fi; 
        end 
    end

end % methods

end % class
