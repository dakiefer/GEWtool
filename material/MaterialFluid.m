classdef MaterialFluid
% MaterialFluid - Represent fluid mechanical material data.
% Subclass of "MaterialIsotropic". It has the same features as the superclass
% but enables to define and load a material with only the density and either the 
% bulk modulus B or the wave velocity cf.
% 
% Example:
% mat = MaterialFluid('water')   % load from water.json (anywhere on path)
% mat = MaterialFluid('anyName', bulkModulus, density);
%
% See also: MaterialIsotropic, Material.
%
% 2025 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France

properties
    name string  % Name of the fluid
    B    double  % Bulk modulus = first Lamé parameter
    rho  double  % Mass density
end

properties (Dependent)
    cl     % Acoustic wave speed
end

methods
     function obj = MaterialFluid(varargin)
        % MATERIALFLUID - Create a fluid material object.
        % 
        % Usage:
        % mat = MaterialFluid('water')   % load from water.json (anywhere on path)
        % mat = MaterialFluid('anyName', bulkModulus, density);
        %
        if nargin == 1 && ischar(varargin{1}) || isstring(varargin{1})  % load by name
            filename = [char(varargin{1}) '.json'];
            data = jsondecode(fileread(filename));
            if ~(isfield(data,'B') && isfield(data,'rho') && isfield(data,'name'))
                error('GEWTOOL:MaterialFluid','Incorrect format of JSON file. The data should contain the fields "name" (descriptive name), "B" (bulk modulus) and "rho" (mass density).');
            end
            varargin = {data.name, data.B, data.rho}; % convert to cell array
        end
        name = varargin{1}; 
        B = varargin{2}; 
        rho = varargin{3}; 
        if nargin >= 4 && string(varargin{4}) == "wavespeed"
            B = MaterialIsotropic.wavespeed2lame(B,0,rho); % B is initially the wave speed
        end 
        obj.B = B; obj.rho = rho; obj.name = string(name); 
     end
    function obj = set.B(obj, B)
        if B <= 0
            error('GEWTOOL:Material:setrange', 'The bulk modulus B needs to be greater than zero.');
        end
        obj.B = B;
    end
    function cl = get.cl(obj)
        cl = sqrt(obj.B/obj.rho);
    end
    function obj = set.cl(obj, cl)
        if cl <= 0
            error('GEWTOOL:Material:setrange', 'The wave speed cl needs to be greater than zero.');
        end
        obj.B = MaterialIsotropic.wavespeed2lame(cl,0,obj.rho);
    end
    function [cs, eu] = wavespeeds(obj, ek)
        if nargin <= 1
            ek = [1;0;0];
        end
        ek = ek(:)/norm(ek);
        cs = zeros(3,1);
        cs(1) = sqrt(obj.B/obj.rho); 
        et = null(ek.'); % svd to compute orthogonal vectors to ek
        eu = [ek, et];
    end

    %% overload operators: 
    function ret = eq(a, b)
        % eq - Test if bulk modulus B and density rho are the same for materials a and b.
        % Usage: 
        % isEq = eq(a, b);
        % isEq = a == b;
        ret = a.B == b.B && a.rho == b.rho;
    end
    function ret = ne(a, b)
        ret = ~eq(a, b);
    end
end

end