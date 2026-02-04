function [u] = displacement(dat)
% DISPLACEMENT returns the displacements u. 
% 
% IMPORTANT: displacement(dat) might be different to dat.u!
% For a piezoelectric plate, displacement(dat) contains only the displacements,
% while dat.u contains all field variables, i.e., it also includes the
% electrostatic potential phi. 
% 
% Usage: 
% > u = displacement(dat);
% 
% Arguments: 
% - dat:    GEWdat object that stores the computed modes. 
%

if ~isscalar(dat) % compute recursively for every waveguide problem in the vector "gew"
    u = arrayfun(@displacement,dat,'UniformOutput',false); % apply to every object in the array "dat"
    return;
end

dof = 1:length(dat.gew.udof); % piezoelectric plate does not include potential in gew.udof

field = dat.u;
u = cell(size(field));
for i = 1:length(u) % for every layer
    u{i} = field{i}(:,:,:,dof);
end

end
