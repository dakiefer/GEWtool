function [cg] = groupVel(gew, dat)
% groupVel - Group velocity magnitude |cg|.
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
% See also: groupVelAxial, energyVelAxial.
% 
% 2022-2024 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France

if isa(gew,"CylinderCircumferential")
    warning('Circumferential waves do not support this function yet. The results might be wrong.');
end
if gew.decouplesLambvsSH
    cg = abs(groupVelAxial(gew, dat));
else
    error('GEWTOOL:groupVel:notimplemented', ...
        'Your waveguide exhibits nonzero transverse group velocity components, which is not implemented yet. You can instead compute energyVelMag() to achieve the exact same result.');
end

end
