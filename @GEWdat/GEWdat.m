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
    cp           % [Nk x Nw] phase velocities in m/s
    ce           % [Nk x Nw] or [Nk x Nw x Nudof] energy velocities in m/s. If Lamb and SH waves decouple, only the axial component of ce is computed. Otherwise, all components are computed.
end % properties

properties (Dependent)
    u        % {Nlay} of [Nk x Nw x N x Nudof] field components (displacements, potentials)
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
        u = fieldComponents(obj);
    end
    function obj = set.u(obj,u)
        obj.Psi = eigenVecs(obj.gew,u);
    end
    function cp = get.cp(obj)
        cp = obj.w./obj.k; 
    end
    function ce = get.ce(obj)
        if obj.gew.decouplesLambvsSH
            ce = energyVelAxial(obj); 
        else
            ce = energyVelVec(obj); 
        end
    end
    function ind = isPropagative(obj,reltol)
        if nargin < 2, reltol = 1e-5; end
        ind = cell(size(obj));
        for i = 1:length(obj)
            krMag = abs(real(obj(i).k*obj(i).gew.h));
            kiMag = abs(imag(obj(i).k*obj(i).gew.h));
            ind{i} = krMag < reltol*kiMag;
        end
    end
    function h = plot(dat, varargin)
        % plot - Plot frequency-wavenumber dispersion curves.
        hasStyle = nargin ~= 1; % if lineSpec not provided, we will be choosing a default
        if ~dat(1).gew.isDissipative && isreal(dat(1).k) % we computed at const k
            wrange = dat(1).w(end,1); 
            style = '-';
        else % we computed at const w
            wrange = dat(1).w(end,end);
            style = '.';
        end
        holdStat = ishold; % creats new axis if not yet existent
        for i = 1:length(dat)
            ww = dat(i).w; kk = real(dat(i).k);
            ww(end+1,:) = nan; kk(end+1,:) = nan; % plote each mode discontinuous
            args = varargin; 
            if ~any(strcmpi(varargin,'DisplayName')) % case insensitive
                args{end+1} = 'DisplayName'; args{end+1} = dat(i).gew.family; 
            end
            if hasStyle
                h(i) = plot(kk(:)/1e3, ww(:)/2/pi/1e6, args{:});
            else
                h(i) = plot(kk(:)/1e3, ww(:)/2/pi/1e6, style, args{:});
            end
            hold on;
        end
        if ~holdStat, hold off; end % restore hold status to previous value
        ylim([0, 1]*wrange/2/pi/1e6); 
        xlabel('wavenumber $k$ in rad/mm','Interpreter','latex'), 
        ylabel('frequency $\omega/2\pi$ in MHz','Interpreter','latex')
        legend('Location', 'southeast');
        if nargout == 0, clear h; end
    end
end

end % class
