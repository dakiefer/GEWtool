function obj = Plate(mats, zs, Ns)
% Plate - Wrapper function to create a PlateClosed or PlateLeaky object
% Arguments: 
% - mats:  materials [1 x Nlay], either of class "Material" or a struct
%          needs to support mats.rho (scalar) and mats.c (3x3x3x3).
% - zs:    coordinates of interfaces in meter [1 x Nlay+1]
%          instead, if zs is of size [1 x Nlay]: zs are interpreted as the 
%          thickness of each layer.
% - Ns:    discretization order for each layer [1 x Nlay]
% 
% Example:
% mat = Material('steel'); % load material data
% h = 1e-3; % thickness in m
% N = 20; % discretization (number of nodal points)
% plate = Plate(mat, h, N); % create waveguide description
% 

% extract exterior half-spaces and crop materials "mats" and coordinates "zs"
if ~iscell(mats)
    mats = num2cell(mats);
end
halfSpaces = PlateLeaky.parseLoading(mats,zs);
mats = mats(~isinf(zs)); zs = zs(~isinf(zs)); % crop to finite layers

% standard specification of materials "mats" and coordinates "zs"
if isscalar(zs) && length(mats) > 1 % use same thickness for all layers
    zs = zs*ones(size(mats)); % expand into a vector of thicknesses
end
if length(mats) == length(zs) % thicknesses have been provided
    if isscalar(zs)
        zs = zs*[-1/2, 1/2]; % single layer has centered coordinate 
    else
        zs = [0, cumsum(zs)]; % coordinates of interfaces starting from 0
    end
elseif length(zs) ~= length(mats)+1
    error('GEWTOOL:Plate:wrongArguments','Provide either a thickness for each layer or the coordinates of the interfaces.');
end

if isempty(halfSpaces)
    obj = PlateClosed(mats, zs, Ns); 
else
    obj = PlateLeaky(mats, zs, Ns, halfSpaces); 
end

end