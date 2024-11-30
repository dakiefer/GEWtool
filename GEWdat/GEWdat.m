classdef GEWdat
% GEWdat - Store and postprocess guided waves solutions.
% Objects of this class are returned by the computeW and computeK functions.
% 
% See also computeW, computeK.
% 
% 2024 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France

properties (Access = public)
    gew          % Waveguide description
	k            % [Nk x Nw] wavenumbers in rad/m
	w            % [Nk x Nw] angular frequencies in rad/s
	Psi = [] 	 % [Nk x Nw x NgdofFree] eigenvectors 
    Nk = 0       % number of wavenumbers (total number of solutions = Nk*Nw) 
    Nw = 0       % number of frequencies (total number of solutions = Nk*Nw) 
    cp           % phase velocity
    ceMag        % energy velocity - magnitude
    cex          % energy velocity - axial component (aligned with k-vector)
    cey          % energy velocity - transverse component (orthogonal to k-vector)
end % properties

properties (Dependent)
    u        % {Nlay} of [Nk x Nw x N x Nudof] unknown field components (displacements, potentials)
end

methods 
    function obj = GEWdat(gew,k,w,Psi)
        % GEWdat - Create GEWdat object. 
        % Arguments: 
        % - gew:   object of class 'Waveguide'
        % - k:     [Nk x Nw] or [Nk x 1] wavenumbers in rad/m
        % - w:     [Nk x Nw] or [1 x Nw] angular frequencies in rad/s
        % - Psi:   [Nk x Nw x NgdofFree] eigenvectors 
        if ~isa(gew,'Waveguide')
            error('GEWtool:GEWdat', 'The first argument needs to be an object of type Waveguide are a subclass thereof.'); 
        end
        obj.gew = copy(gew);
        obj.k = k; 
        obj.Nk = size(k,1);
        obj.w = w; 
        obj.Nw = size(w,2);
        if nargin > 3
            obj.Psi = Psi;
        end
    end
    function k = get.k(obj)
        % get.k - getter function for obj.k. 
        % If k is a column vector, it will be expanded to the same size as obj.w.
        if iscolumn(obj.k) && obj.Nw > 1
            k = obj.k.*ones([obj.Nk, obj.Nw]);
        else 
            k = obj.k; 
        end
    end
    function w = get.w(obj)
        % get.w - getter function for obj.w. 
        % If w is a row vector, it will be expanded to the same size as obj.k.
        if isrow(obj.w) && obj.Nk > 1
            w = obj.w.*ones([obj.Nk, obj.Nw]);
        else 
            w = obj.w;
        end
    end
    function u = get.u(obj)
        u = unknowns(obj);
    end
    function obj = set.u(obj,u)
        obj.Psi = eigenVecs(obj.gew,u);
    end
    function cp = get.cp(obj)
        cp = obj.w./obj.k; 
    end
    function ceMag = get.ceMag(obj)
        ceMag = energyVelMag(obj); 
    end
    function cex = get.cex(obj)
        cex = energyVelAxial(obj); 
    end
    function cey = get.cey(obj)
        cey = energyVelTransverse(obj); 
    end
    % u = unknowns(obj)
end

end % class
