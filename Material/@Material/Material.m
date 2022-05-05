classdef Material
% MATERIAL represent and load material data.
%
% 2022 - Daniel Kiefer
% Institut Langevin, Paris, France

properties
    name        % name for the material: string
    symmetry    % symmetry class: e.g., "isotropic": string
    c double    % 4th order stiffness tensor
    rho double  % mass density [1x1]
end

properties (Dependent)
    C double % Voigt notated stiffness [6x6]
    cl       % (quasi)-longitudinal waves speed
    ct       % (quasi)-transverse waves speed 1
    ct2      % (quasi)-transverse waves speed 2
end

methods
    function obj = Material(varargin)
        if nargin == 1 && (isa(varargin{1}, 'Material') || isa(varargin{1}, 'struct')) % convert from subclasses to superclass
            data = varargin{1};
        elseif nargin == 1 && (ischar(varargin{1}) || isstring(varargin{1}))  % load by name
            matname = varargin{1};
            filename = [matname '.json'];
            data = jsondecode(fileread(filename));
            if ~isfield(data, 'symmetry') % add field if not in json file
                if isfield(data, 'lambda') && isfield(data, 'mu')
                    data.symmetry = 'isotropic';
                else
                    data.symmetry = 'unknown'; % might still be isotropic... 
                end
            end
        elseif nargin == 3  % name, C or c, rho 
            validateattributes(varargin{1},{'char'},{'vector'},1); % name
            data.name = varargin{1};
            if ~isnumeric(varargin{2}) % stiffness
                error('Expected stiffness (argument 2) to be numeric.');
            elseif isequal(size(varargin{2}), [6,6]) % Voigt noteted tensor
                data.C = varargin{2};
            elseif isequal(size(varargin{2}), [3,3,3,3]) % 4th order tensor
                data.c = varargin{2};
            else
                error('Expected argument 2 to be of size [6x6] or [3x3x3x3]');
            end
            validateattributes(varargin{3},{'numeric'},{'size',[1 1]},3); % density
            data.rho = varargin{3};
            data.symmetry = 'anisotropic';
        elseif nargin == 4 && (ischar(varargin{1}) || isstring(varargin{1})) % name, lambda, mu and rho 
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
            data.c = data.lambda*II + data.mu*(permute(II, [1 3 4 2]) + permute(II, [1 3 2 4]));
        elseif isfield(data, 'C') 
            data.c = voigt2tensor(data.C);
        end
        obj.c = data.c;
        obj.rho = data.rho;
        obj.name = data.name;
        obj.symmetry = data.symmetry;
    end

    function C = get.C(obj)
        C = tensor2voigt(obj.c);
    end
    function obj = set.C(obj, C)
        obj.c = voigt2tensor(C);
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
        ret = ~any(a.c ~= b.c, 'all') && a.rho == b.rho;
    end
    function ret = ne(a, b)
        ret = ~eq(a, b);
    end
    function plot(varargin)
        plotSlownessCurve(varargin{:});
    end
    
end % methods
end % class
