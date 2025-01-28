function [cg] = groupVelAxial(dat)
% groupVelAxial - Group velocity component cg_x = dw/dk along the wave vector k.
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
% Usage: 
% > cgx = groupVelAxial(dat);
% 
% Arguments: 
% - dat:    GEWdat object that stores the computed modes. 
%
% See also: energyVel, energyVelVec.
% 
% 2022-2024 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France

if ~isscalar(dat) % compute recursively for every waveguide problem in the vector "dat"
    cg = arrayfun(@groupVelAxial,dat,'UniformOutput',false); % apply to every object in the arrays "dat"
    return;
end

gew = dat.gew;
if gew.isDissipative
    warning('GEWTOOL:groupVel:dissipative', ...
        'The waveguide is dissipative and the group velocity is not meaningful. You should compute the energy velocity instead. I will proceed anyways.');
end
if any(~isreal(dat.k))
    warning('GEWTOOL:groupVel:complex', ...
        'The group velocity of modes with complex wavenumbers is meaningless, use energyVelAxial() instead. I will proceed anyways and real modes will yield a meaningful group velocity.');
end
M = gew.op.M; % mass matrix 
L2 = gew.op.L2; % stiffness matrix in k^2
L1 = gew.op.L1; % stiffness matrix in k
w = dat.w/gew.np.fh0*gew.np.h0; % normalize frequency like the above operators
k = dat.k*gew.np.h0; % wavenumbers
Psi = dat.Psi;

cg = zeros(size(k)); 
for i = 1:size(k,1)
    for j = 1:size(k,2)
        Psi0 = squeeze(Psi(i,j,:));
        cg(i, j) = (Psi0'*(2*k(i,j).*L2 - 1i*L1)*Psi0) ./ (2*w(i,j)*Psi0'*M*Psi0);
    end
end
cg = cg*gew.np.fh0;

end
