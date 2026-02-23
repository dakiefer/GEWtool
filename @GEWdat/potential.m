function [phi] = potential(dat)
% POTENTIAL returns the electrostatic potential phi. 
% 
% Only applicable to piezoelectric waveguides. A zero-cell array of appropriate size 
% is returned otherwise. 
% 
% Usage: 
% > phi = potential(dat);
% 
% Arguments: 
% - dat:    GEWdat object that stores the computed modes. 
%

if ~isscalar(dat) % compute recursively for every waveguide problem in the vector "gew"
    phi = arrayfun(@potential,dat,'UniformOutput',false); % apply to every object in the array "dat"
    return;
end

field = dat.u;
phi = cell(size(field));

if isa(dat.gew.lay{1}, 'LayerPlatePiezo') % extract potentials
    dofPhi = length(dat.gew.udof)+1; % potentials are stored after all displacement components 
                                 % piezoelectric plate does not include potential in gew.udof
    for i = 1:length(phi) % for every layer
        phi{i} = field{i}(:,:,:,dofPhi);
    end
else % no potentials -> set to zero
    s = size(field{1});
    for i = 1:length(phi)
        phi{i} = zeros(s(1:3)); % initialize potential to zero for layers without piezoelectric properties
    end
end

end
