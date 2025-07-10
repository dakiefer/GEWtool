function obj = Plate(mats, zs, Ns)
% Plate - Wrapper function to create a PlateClosed or PlateLeaky object
% Arguments: 
% - mats:  materials [1 x Nlay], either of class "Material" or a struct
%          needs to support mats.rho (scalar) and mats.c (3x3x3x3).
% - zs:    thickness of each layer in meter [1 x Nlay] or,
%          coordinates of interfaces in meter [1 x Nlay+1]
%          -> infinite thickness or coordinate is interpreted as a half-space.
% - Ns:    discretization order for each layer [1 x Nlay]
% 
% Example one layer:
% mat = Material('steel'); % load material data
% h = 1e-3; % thickness in m
% N = 20; % discretization (number of nodal points)
% plate = Plate(mat, h, N); % create waveguide description
% 
% Example two layers:
% mat1 = Material('steel'); % load material data
% mat2 = Material('pmma'); 
% h = 1e-3; % thickness in m
% N = 20; % discretization (number of nodal points)
% plate = Plate({mat1 mat2}, h, N); % two layers of same thickness h
%
% Example leaky plate:
% steel = Material('steel'); % load material data
% water = Material('water'); 
% h = 1e-3; % thickness in m
% N = 20; % discretization (number of nodal points)
% plate = Plate({steel water}, [h inf], N); % steel plate with water on top
% 
% 2025 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France

    % convert material list to a cell array 
    if ~iscell(mats)
        mats = num2cell(mats);
    end
    
    % convert to list of thicknesses "hs" for each layer (if not yet the case):
    if isscalar(zs) && length(mats) > 1 % one thickness provided for multiple layers
        zs = zs*ones(size(mats)); % use same thickness for all layers
    end
    if length(zs) == length(mats) + 1 % coordinates of interfaces provided
        zs = diff(zs); % convert to thickness for each layer
    end
    hs = zs; % here we will always have a list of thicknesses for each mats{i}

    % extract half-spaces and crop materials "mats" and coordinates "zs"
    halfSpaces = PlateLeaky.parseLoading(mats,hs);
    mats = mats(~isinf(hs)); hs = hs(~isinf(hs)); % crop to finite layers
    
    % convert to coordinates "zs" for all layers of finite thickness:
    if isscalar(hs)
        zs = hs*[-1/2, 1/2]; % single layer has centered coordinate 
    else
        zs = [0, cumsum(hs)]; % coordinates of interfaces starting from 0
    end
    
    % error checking:
    if length(zs) ~= length(mats)+1
        error('GEWTOOL:Plate:wrongArguments','Provide either a thickness for each layer or the coordinates of the interfaces.');
    end
    
    % create "PlateClosed" or "PlateLeaky" object:
    if isempty(halfSpaces) % no exterior half-spaces
        obj = PlateClosed(mats, zs, Ns); % closed plate (no radiation)
    else
        obj = PlateLeaky(mats, zs, Ns, halfSpaces); 
    end

end
