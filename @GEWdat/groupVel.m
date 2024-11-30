function [cg] = groupVel(dat)
% groupVel - Group velocity real-part magnitude |Re{cg}|.
% 
% The axial group velocity is computed based on the eigenvectors. No explicit
% differntiation of the solutions is performed, i.e., the accuracy is
% independent of the discretization of your k-w axes and the computation is
% performed for every single mode independently. 
% 
% The group velocity is only meaningul for nondissipative waveguides. In a more
% general setting, you can compute the energy velocity, which for the
% nondissipative case is identical to the group velocity.
% 
% The imaginary part is considered meaningless. After taking the magnitude it
% can no longer be separated from the real part, so we remove it before. 
% 
% See also: groupVelAxial, energyVelAxial.
% 
% 2022-2024 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France

if ~isscalar(dat) % compute recursively for every waveguide problem in the vector "dat"
    cg = arrayfun(@groupVel,dat,'UniformOutput',false); % apply to every object in the arrays "dat"
    return; 
end

if dat.gew.decouplesLambvsSH
    cg = abs(real(groupVelAxial(dat)));
else
    error('GEWTOOL:groupVel:notimplemented', ...
        'Your waveguide exhibits nonzero transverse group velocity components, which is not implemented yet. You can instead compute energyVelMag() to achieve the exact same result or you can compute the axial component of the group velocity vector.');
end

end
