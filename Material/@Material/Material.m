classdef Material
% MATERIAL represent and load material data.
%
% 2022 - Daniel Kiefer
% Institut Langevin, Paris, France

properties
    name        % name for the material: string
    symmetry    % symmetry class: e.g., "isotropic": string
    C double    % Voigt notated stiffness [6x6]
    rho double  % mass density [1x1]
end

properties (Dependent)
    c double % 4th order stiffness tensor
    cl       % (quasi)-longitudinal waves speed
    ct       % (quasi)-transverse waves speed 1
    ct2      % (quasi)-transverse waves speed 2
end

methods
    function obj = Material(varargin)
        if nargin == 1 && (ischar(varargin{1}) || isstring(varargin{1}))  % load by name
            matname = varargin{1};
            [dir, ~] = fileparts(which('Material'));
            dir = fileparts(dir); % remove the @Material folder
            dir = fullfile(dir, 'database');
            filename = fullfile(dir, [matname '.json']);
            data = jsondecode(fileread(filename));
        elseif nargin == 3 && (ischar(varargin{1}) || isstring(varargin{1})) 
            validateattributes(varargin{2},{'numeric'},{'size',[6 6]});
            validateattributes(varargin{3},{'numeric'},{'size',[1 1]});
            data.name = varargin{1};
            data.C = varargin{2};
            data.rho = varargin{3};
            data.symmetry = 'anisotropic';
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

    function c = get.c(obj)
        c = voigt2tensor(obj.C);
    end
    function obj = set.c(obj, c)
        obj.C = tensor2voigt(c);
    end
    function cl = get.cl(obj)
        cs = wavespeeds(obj);
        cl = cs(1);
    end
    function ct = get.ct(obj)
        cs = wavespeeds(obj);
        ct = cs(2);
    end
    function ct2 = get.ct2(obj)
        cs = wavespeeds(obj);
        ct2 = cs(3);
    end

    function obj = permute(obj, perm)
        % PERMUTE Permute the material coordinate system. The default
        % transformation is ex-ey-ez -> ez-ex-ey, i.e. 
        % 3->1, 1->2, 2->3, 6->4, 4->5, 5->6, given by perm = [3,1,2,6,4,5].
        if nargin < 2
            perm = [3,1,2,6,4,5]; % 3->1, 1->2, 2->3, 6->4, 4->5, 5->6
        end
        obj.C = obj.C(perm,perm);
    end

    plotSlownessCurve(varargin)
    [cs, eu] = wavespeeds(obj, ek)

    %% overload operators: 
    function ret = eq(a, b)
        ret = ~any(a.C ~= b.C, 'all') && a.rho == b.rho;
    end
    function ret = ne(a, b)
        ret = ~eq(a, b);
    end
    function plot(varargin)
        plotSlownessCurve(varargin{:});
    end
    
end % methods
end % class
