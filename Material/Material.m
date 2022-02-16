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
end % methods
end % class