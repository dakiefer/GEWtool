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
% 2022 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France

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
end % methods

end % class
