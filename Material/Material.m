classdef Material < matlab.mixin.Copyable
% MATERIAL load and represent material data.
%
% 2022 - Daniel Kiefer
% Institut Langevin, Paris, France

properties
    name
    symmetry
    c double
    rho double
end

properties (Dependent)
    C double 
    cl
    ct
end

methods
    function obj = Material(varargin)
        if nargin == 1 && (ischar(varargin{1}) || isstring(varargin{1}))  % load by name
            matname = varargin{1};
            [dir, ~] = fileparts(which('Material'));
            dir = fullfile(dir, 'database');
            filename = fullfile(dir, [matname '.json']);
            data = jsondecode(fileread(filename));
        elseif nargin == 4 && (ischar(varargin{1}) || isstring(varargin{1})) 
            validateattributes(varargin{2},{'numeric'},{'size',[1 1]});
            validateattributes(varargin{3},{'numeric'},{'size',[1 1]});
            validateattributes(varargin{4},{'numeric'},{'size',[1 1]});
            data.name = varargin{1}; 
            data.lambda = varargin{2};
            data.mu = varargin{3};
            data.rho = varargin{4};
            data.symmetry = 'isotropic';
        end
        if strcmp(data.symmetry, 'isotropic') 
            II = eye(3).*shiftdim(eye(3), -2); % 4th order "unit tensor"
            obj.c = data.lambda*II + data.mu*(permute(II, [1 3 4 2]) + permute(II, [1 3 2 4]));
        else 
            obj.c = voigt2tensor(data.C);
        end
        obj.rho = data.rho;
        obj.name = data.name;
        obj.symmetry = data.symmetry;
    end
    
    function [cs, eu] = wavespeeds(obj, ek)
        % WAVESPEEDS Compute the wave speeds by solving the Kelvin-Christoffel
        % equation (eigenvalue problem for cl, ct1, ct2).
        if nargin < 2, ek = [1; 0; 0]; end
        ek = ek(:); ekS = shiftdim(ek, -3);
        D = 1/obj.rho*sum(sum(ek.*obj.c.*ekS, 4), 1); % Kristoffel tensor: 1/rho ek.c.ek
        D = squeeze(D); % 3x3 matrix
        if nargout == 1
            cs = sort(sqrt(eig(D)), 'descend');
        else
            [eu, cs2] = eig(D, 'vector'); 
            [cs, ind] = sort(sqrt(cs2), 'descend');
            eu = eu(:,ind);
        end
    end
    
    function cl = get.cl(obj)
        cs = wavespeeds(obj);
        cl = cs(1);
    end
    function ct = get.ct(obj)
        cs = wavespeeds(obj);
        ct = cs(2:3);
    end
    
end % methods
end % class